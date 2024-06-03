import argparse
import os
import sys
import numpy as np
import warnings
import torch
import logging
import pandas as pd

from model import TFBSNet
from utils import set_seed

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
    parser.add_argument("--data", type=str, default=None)
    parser.add_argument("--model", type=str, default=None)
    parser.add_argument("--assay", type=str, default='atac')
    parser.add_argument("--out_dir", type=str,
                        default=None, help="Output directory")
    parser.add_argument("--out_name", type=str,
                        default=None, help="Output name")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    return parser.parse_args()


def main():
    args = parse_args()

    set_seed(args.seed)

    logging.info("Loading input files")
    data = np.load(args.data)

    # Setup model
    logging.info("Creating model")
    if args.assay == 'atac':
        model = TFBSNet(n_channels=6)
    elif args.assay == 'access_atac':
        model = TFBSNet(n_channels=8)

    state_dict = torch.load(args.model)
    model.load_state_dict(state_dict["state_dict"])
    device = torch.device("cuda")
    model.to(device)
    model.eval()

    logging.info("Predicting started")
    preds = np.empty(shape=data['x_test'].shape[0])
    for i, x in enumerate(data['x_test']):
        x = torch.tensor(x)
        x = x.unsqueeze(dim=0)
        pred = model(x.to(device)).detach().cpu().view(-1)

        # convert to probability
        preds[i] = torch.sigmoid(pred).numpy()[0]

    logging.info(f"Predicting finished")

    df = pd.DataFrame(data={"target": data['y_test'],
                            "pred": preds})

    df.to_csv(f'{args.out_dir}/{args.out_name}.csv', index=False)


if __name__ == "__main__":
    main()
