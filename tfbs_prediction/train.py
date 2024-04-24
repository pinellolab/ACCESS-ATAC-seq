import argparse
import os
import sys
import numpy as np
import warnings
import torch
import logging
from torch.optim import Adam
from torch.utils.tensorboard import SummaryWriter
from torch.optim.lr_scheduler import ReduceLROnPlateau

from model import TFBSNet
from utils import set_seed
from dataset import get_dataloader

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)

if not sys.warnoptions:
    warnings.simplefilter("ignore")


def parse_args():
    parser = argparse.ArgumentParser()

    # Required parameters
    parser.add_argument("--train_data", type=str, default=None)
    parser.add_argument("--valid_data", type=str, default=None)
    parser.add_argument("--test_data", type=str, default=None)
    parser.add_argument("--assay", type=str, default='atac')
    parser.add_argument(
        "--epochs", type=int, default=200, help="Number of epochs for training"
    )
    parser.add_argument("--log_dir", type=str, default=None, help="Log directory")
    parser.add_argument("--log_name", type=str, default=None, help="Log name")
    parser.add_argument("--out_dir", type=str, default=None, help="Output directory")
    parser.add_argument("--out_name", type=str, default=None, help="Output name")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--batch_size", type=int, default=48, help="Batch size")
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
    train_data = np.load(args.train_data)
    valid_data = np.load(args.valid_data)
    test_data = np.load(args.test_data)
    
    train_dataloader = get_dataloader(
        x=train_data['x'],
        y=train_data['y'],
        batch_size=args.batch_size,
        drop_last=True,
        shuffle=True,
        train=True,
    )
    valid_dataloader = get_dataloader(
        x=valid_data['x'],
        y=valid_data['y'],
        batch_size=args.batch_size,
        drop_last=False,
        shuffle=False,
        train=True,
    )

    test_dataloader = get_dataloader(
        x=test_data['x'],
        y=test_data['y'],
        batch_size=args.batch_size,
        drop_last=False,
        shuffle=False,
        train=True,
    )    

    # Setup model
    logging.info("Creating model")
    if args.assay == 'atac':
        model = TFBSNet(n_channels=5)
    elif args.assay == 'access_atac':
        model = TFBSNet(n_channels=6)
        
    device = torch.device("cuda")
    model.to(device)

    # Setup loss and optimizer
    criterion = torch.nn.CrossEntropyLoss()
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
