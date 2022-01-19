"""Utility functions for python scripts."""
import argparse


def print_args(args: argparse.Namespace, fmt: str = ">40"):
    """Print args in a readable format."""
    for arg, value in vars(args).items():
        try:
            print(f"{arg:{fmt}}: {value:{fmt}}")
        except Exception:
            print(f"{arg:{fmt}}: {'None (NoneType)':{fmt}}")
