import argparse
import os
import sys
import numpy as np
import pandas as pd
import warnings
import torch
import logging
from torch.optim import Adam
import pyBigWig
import pyranges as pr
import pysam

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

chrom_valid = ["chr5", "chr8", "chr20"]
chrom_train = [
    "chr2",
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
        "--peak_file", type=str, default=None, help="BED file of input genomic regions"
    )
    parser.add_argument(
        "--ref_fasta", type=str, default=None, help="FASTQ file for reference genome"
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="random seed for initialization"
    )
    parser.add_argument("--batch_size", type=int, default=128, help="Batch size")
    return parser.parse_args()


def train(model, dataloader, criterion, optimizer, device):
    model.train()

    train_loss = 0.0
    for x, y in dataloader:
        pred = model(x.to(device))
        loss = criterion(pred.view(-1), y.to(device))

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
        loss = criterion(pred.view(-1), y.to(device))
        valid_loss += loss.item() / len(dataloader)

    return valid_loss


def main():
    args = parse_args()

    set_seed(args.seed)

    logging.info("Loading input files")
    bw = pyBigWig.open(args.bw_file)
    grs = pr.read_bed(args.peak_file)
    fasta = pysam.FastaFile(args.ref_fasta)

    logging.info("Preprocessing input regions")
    grs = pad_and_split(grs, k=128)

    # split the regions for training and validation
    grs_train = grs[grs.Chromosome.isin(chrom_train)]
    grs_valid = grs[grs.Chromosome.isin(chrom_valid)]

    logging.info("Generating data for training and validation")
    train_x = np.empty(shape=(len(grs_train), 128, 4))
    train_y = np.empty(shape=(len(grs_train), 128))
    for i, (chrom, start, end) in enumerate(
        zip(grs_train.Chromosome, grs_train.Start, grs_train.End)
    ):
        # get DNA sequence
        seq = str(fasta.fetch(chrom, start, end)).upper()

        # convert DNA sequence to one-hot encode
        train_x[i] = one_hot_encode(seq=seq)
        train_y[i] = np.array(bw.values(chrom, start, end))

    valid_x = np.empty(shape=(len(grs_valid), 128, 4))
    valid_y = np.empty(shape=(len(grs_valid), 128))
    for i, (chrom, start, end) in enumerate(
        zip(grs_valid.Chromosome, grs_valid.Start, grs_valid.End)
    ):
        # get DNA sequence
        seq = str(fasta.fetch(chrom, start, end)).upper()

        # convert DNA sequence to one-hot encode
        valid_x[i] = one_hot_encode(seq=seq)
        valid_y[i] = np.array(bw.values(chrom, start, end))

    # logging(f"Training data: {train_x.shape[0]}, validation data: {valid_x.shape[0]}")

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
    model = BiasNet(seq_len=128)
    device = torch.device("cuda")
    model.to(device)

    # Setup loss and optimizer
    criterion = torch.nn.MSELoss()
    optimizer = Adam(model.parameters(), lr=3e-04, weight_decay=1e-4)

    logging.info("Training started")
    for epoch in range(10):
        train_loss = train(
            dataloader=train_dataloader,
            model=model,
            criterion=criterion,
            optimizer=optimizer,
            device=device,
        )
        
        print(train_loss)


if __name__ == "__main__":
    main()
