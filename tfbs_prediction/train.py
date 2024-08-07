import argparse
import os
import sys
import numpy as np
import pandas as pd
import warnings
import torch
import logging
from torch.optim import Adam
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
    parser.add_argument("--assay", type=str, default='atac')
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--model_path", type=str, default=None)
    parser.add_argument("--log_path", type=str, default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--batch_size", type=int, default=48)
    return parser.parse_args()


def train(model, dataloader, criterion, optimizer, device):
    model.train()

    train_loss = 0.0
    for x, y in dataloader:
        pred = model(x.to(device)).view(-1)
        loss = criterion(pred.float(), y.to(device).float())

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        train_loss += loss.item() / len(dataloader)

    return train_loss


def valid(model, dataloader, criterion, device):
    model.eval()

    valid_loss = 0.0
    for x, y in dataloader:
        pred = model(x.to(device)).view(-1)
        loss = criterion(pred.float(), y.to(device).float())

        valid_loss += loss.item() / len(dataloader)

    return valid_loss


def main():
    args = parse_args()

    set_seed(args.seed)

    logging.info("Loading input files")
    train_data = np.load(args.train_data)
    valid_data = np.load(args.valid_data)

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

    # Setup model
    logging.info(f"Creating model for {args.assay}")
    if args.assay == 'atac' or args.assay == 'access':
        model = TFBSNet(n_channels=6)
    elif args.assay == 'access_atac':
        model = TFBSNet(n_channels=8)
    elif args.assay == 'dna':
        model = TFBSNet(n_channels=4)

    device = torch.device("cuda")
    model.to(device)

    # Setup loss and optimizer
    criterion = torch.nn.BCEWithLogitsLoss()
    optimizer = Adam(model.parameters(), lr=1e-04, weight_decay=1e-4)
    scheduler = ReduceLROnPlateau(optimizer, "min", min_lr=1e-5, patience=5)

    """ Train the model """
    logging.info("Training started")
    best_score = np.inf

    epochs, train_losses, valid_losses, best_scores = [], [], [], []
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

        # save model if find a better validation score
        if valid_loss < best_score:
            best_score = valid_loss
            state = {
                "state_dict": model.state_dict(),
                "train_loss": train_loss,
                "valid_loss": valid_loss,
                "epoch": epoch,
            }
            torch.save(state, args.model_path)
            # Reset patience counter
            patience = 10
        else:
            # early stop
            patience -= 1
            if patience == 0:
                logging.info("Early stop!")
                break

        logging.info(
            f"epoch: {epoch}, train: {train_loss:.5f}, valid: {valid_loss:.5f}, best: {best_score:.5f}")
        scheduler.step(valid_loss)

        epochs.append(epoch)
        train_losses.append(train_loss)
        valid_losses.append(valid_loss)
        best_scores.append(best_score)

    df = pd.DataFrame(data={"epoch": epochs,
                            "train_loss": train_losses,
                            "valid_loss": valid_losses,
                            "best_loss": best_scores})

    df.to_csv(args.log_path, index=False)
    logging.info(f"Training finished")


if __name__ == "__main__":
    main()
