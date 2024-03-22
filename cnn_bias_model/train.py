import argparse
import os
import sys
import numpy as np
import warnings
import torch
import logging
from torch.optim import Adam
from torch.nn.utils.clip_grad import clip_grad_value_
import pyBigWig
import pyranges as pr
import pysam
from torch.utils.tensorboard import SummaryWriter
from torch.optim.lr_scheduler import ReduceLROnPlateau, StepLR

from model import BiasNet
from utils import pad_and_split, one_hot_encode, set_seed
from dataset import get_dataloader

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)

if not sys.warnoptions:
    warnings.simplefilter("ignore")

chrom_valid = ["chr2", "chr20"]
chrom_train = [
    "chr1",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr21",
    "chr22",
    "chrX",
]


def parse_args():
    parser = argparse.ArgumentParser()

    # Required parameters
    parser.add_argument(
        "--bw_file",
        type=str,
        default=None,
        help="Bigwig file containing the raw signal for training",
    )
    parser.add_argument(
        "--train_regions",
        type=str,
        default=None,
        help="BED file containing genomic regions for training",
    )
    parser.add_argument(
        "--valid_regions",
        type=str,
        default=None,
        help="BED file containing genomic regions for validation",
    )
    parser.add_argument(
        "--ref_fasta", type=str, default=None, help="FASTQ file for reference genome"
    )
    parser.add_argument("--k", type=int, default=128, help="Input sequence size")
    parser.add_argument(
        "--epochs", type=int, default=200, help="Number of epochs for training"
    )
    parser.add_argument("--log_dir", type=str, default=None, help="Log directory")
    parser.add_argument("--log_name", type=str, default=None, help="Log name")
    parser.add_argument("--out_dir", type=str, default=None, help="Output directory")
    parser.add_argument("--out_name", type=str, default=None, help="Output name")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--batch_size", type=int, default=512, help="Batch size")
    return parser.parse_args()


def train(model, dataloader, criterion, optimizer, device):
    model.train()

    train_loss = 0.0
    for x, y in dataloader:
        pred = model(x.to(device))
        loss = criterion(pred, y.to(device))

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        train_loss += loss.item() / len(dataloader)

    return train_loss


def valid(model, dataloader, criterion, device):
    model.eval()

    valid_loss = 0.0
    for x, y in dataloader:
        pred = model(x.to(device))
        loss = criterion(pred, y.to(device))
        valid_loss += loss.item() / len(dataloader)

    return valid_loss


def main():
    args = parse_args()

    set_seed(args.seed)

    logging.info("Loading input files")
    bw = pyBigWig.open(args.bw_file)
    fasta = pysam.FastaFile(args.ref_fasta)

    grs_train = pr.read_bed(args.train_regions)
    grs_valid = pr.read_bed(args.valid_regions)

    logging.info("Preprocessing input regions")
    grs_train = pad_and_split(grs_train, fasta_file=args.ref_fasta, k=args.k)
    grs_valid = pad_and_split(grs_valid, fasta_file=args.ref_fasta, k=args.k)

    # split the regions for training and validation
    # grs_train = grs[grs.Chromosome.isin(chrom_train)]
    # grs_valid = grs[grs.Chromosome.isin(chrom_valid)]

    logging.info("Generating data for training and validation")
    train_x = np.empty(shape=(len(grs_train), args.k, 4), dtype=np.float32)
    train_y = np.empty(shape=(len(grs_train), args.k), dtype=np.float32)
    for i, (chrom, start, end) in enumerate(
        zip(grs_train.Chromosome, grs_train.Start, grs_train.End)
    ):
        # get DNA sequence
        if start < 0:
            continue
        
        seq = str(fasta.fetch(chrom, start, end)).upper()

        # convert DNA sequence to one-hot encode
        train_x[i] = one_hot_encode(seq=seq)
        signal = np.array(bw.values(chrom, start, end))
        signal[np.isnan(signal)] = 0
        train_y[i] = signal

    valid_x = np.empty(shape=(len(grs_valid), args.k, 4), dtype=np.float32)
    valid_y = np.empty(shape=(len(grs_valid), args.k), dtype=np.float32)
    for i, (chrom, start, end) in enumerate(
        zip(grs_valid.Chromosome, grs_valid.Start, grs_valid.End)
    ):
        # get DNA sequence
        seq = str(fasta.fetch(chrom, start, end)).upper()

        # convert DNA sequence to one-hot encode
        valid_x[i] = one_hot_encode(seq=seq)
        signal = np.array(bw.values(chrom, start, end))
        signal[np.isnan(signal)] = 0
        valid_y[i] = signal

    logging.info(f"Number of training: {len(grs_train)}")
    logging.info(f"Number of validation: {len(grs_valid)}")

    train_dataloader = get_dataloader(
        x=train_x,
        y=train_y,
        batch_size=args.batch_size,
        drop_last=True,
        shuffle=True,
        train=True,
    )
    valid_dataloader = get_dataloader(
        x=valid_x,
        y=valid_y,
        batch_size=args.batch_size,
        drop_last=False,
        shuffle=False,
        train=True,
    )

    # Setup model
    logging.info("Creating model")
    model = BiasNet(seq_len=args.k)
    device = torch.device("cuda")
    model.to(device)

    # Setup loss and optimizer
    criterion = torch.nn.MSELoss()
    optimizer = Adam(model.parameters(), lr=1e-04, weight_decay=1e-4)
    scheduler = ReduceLROnPlateau(optimizer, "min", min_lr=1e-5, patience=10)

    """ Train the model """
    os.makedirs(args.log_dir, exist_ok=True)
    log_dir = os.path.join(args.log_dir, f"{args.log_name}")
    tb_writer = SummaryWriter(log_dir=log_dir)

    logging.info("Training started")
    model_path = os.path.join(args.out_dir, f"{args.out_name}.pth")
    best_score = np.Inf
    for epoch in range(args.epochs):
        train_loss = train(
            dataloader=train_dataloader,
            model=model,
            criterion=criterion,
            optimizer=optimizer,
            device=device,
        )
        valid_loss = valid(
            dataloader=valid_dataloader, model=model, criterion=criterion, device=device
        )

        # save log
        tb_writer.add_scalar("Training loss", train_loss, epoch)
        tb_writer.add_scalar("Valid loss", valid_loss, epoch)
        tb_writer.add_scalar("Learning rate", optimizer.param_groups[0]["lr"], epoch)

        # save model if find a better validation score
        if valid_loss < best_score:
            best_score = valid_loss
            logging.info(f"epoch: {epoch}, best score: {best_score}")
            state = {
                "state_dict": model.state_dict(),
                "train_loss": train_loss,
                "valid_loss": valid_loss,
                "epoch": epoch,
            }
            torch.save(state, model_path)
        scheduler.step(valid_loss)

    logging.info(f"Training finished")


if __name__ == "__main__":
    main()
