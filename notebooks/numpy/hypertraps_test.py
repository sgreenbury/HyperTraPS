"""Tests for HyperTraPS class."""
import numpy as np
from hypertraps import HyperTraPS
import pytest


def make_test_data_1(length=7, repeats=1):
    """Left to right binary data."""
    start = np.zeros((length, length))
    end = start.copy()
    for i in range(0, length):
        for j in range(0, i + 1):
            end[i, j] = 1
    start = np.repeat(start, repeats, axis=0).astype(int)
    end = np.repeat(end, repeats, axis=0).astype(int)
    return np.stack([start, end], axis=1)


def make_test_data_2(length=7, repeats=1):
    """Left to right and right to left binary data."""
    X = make_test_data_1(length, repeats=repeats)
    start1, end1 = X[:, 0, :], X[:, 1, :]
    start2 = start1[:, ::-1]
    end2 = end1[:, ::-1]
    start = np.concatenate([start1, start2], axis=0)
    end = np.concatenate([end1, end2], axis=0)
    return np.stack([start, end], axis=1)


def make_zeros_logits(length):
    return np.zeros((length, length))


def make_left_right_logits(length, width=20):
    logits = np.zeros((length, length))
    logits[0, 0] = width
    for i in range(1, length):
        logits[i - 1, i] = width
    return logits


def make_left_right_right_left_logits(length, width=40):
    logits = np.zeros((length, length))
    logits[0, 0] = width // 2
    for i in range(1, length):
        logits[i - 1, i] = width
    logits[-1, -1] = width // 2
    for i in range(1, length):
        logits[i, i - 1] = width
    return logits


# Some lengths of sequences over which we can test the likelihood
# TODO: add test cases and handling for negative and large length
lengths = [1, 4, 5, 7, 10, 20, 30, 40, 47]


@pytest.mark.parametrize(
    ("length", "logits_type"),
    list(
        zip(2 * lengths, len(lengths) * ["zeros"] + len(lengths) * ["optimal"])
    ),
)
def test_left_right_zeros(length, logits_type):
    """Test left to right and right to acquiring binary digits with zeros."""
    # Initialize the left-right/right-left data
    X = make_test_data_1(length=length, repeats=1)

    if logits_type == "zeros":
        # Make zero logits
        logits = np.zeros((length, length))
    elif logits_type == "optimal":
        # Make zero logits
        logits = make_left_right_logits(length)

    # Make HyperTraPS object
    hypertraps = HyperTraPS(logits=logits, verbose=False, n_samples=200)
    actual = hypertraps.log_prob(X)

    # E.g. for length = 7
    # np.log(1 / 7),
    # np.log(2 / 7 * 1 / 6),
    # np.log(3 / 7 * 2 / 6 * 1 / 5),
    # np.log(4 / 7 * 3 / 6 * 2 / 5 * 1 / 4),
    # np.log(5 / 7 * 4 / 6 * 3 / 5 * 2 / 4 * 1 / 3),
    # np.log(6 / 7 * 5 / 6 * 4 / 5 * 3 / 4 * 2 / 3 * 1 / 2),

    # Expected log likelihoods
    if logits_type == "zeros":
        expected = np.array(
            [
                np.sum(np.log(np.arange(1, i + 1)))
                - np.sum(np.log(np.arange(length, length - i, -1)))
                for i in range(1, length + 1)
            ]
        )[:, np.newaxis]
        assert np.allclose(actual, expected)
    elif logits_type == "optimal":
        expected = 0.0
        print(actual, expected)
        assert np.isclose(np.sum(actual), expected, atol=1e-3)


@pytest.mark.parametrize(
    ("length", "logits_type"),
    list(
        zip(2 * lengths, len(lengths) * ["zeros"] + len(lengths) * ["optimal"])
    ),
)
def test_left_right_right_left_zeros(length, logits_type):
    """Test left to right and right to acquiring binary digits with zeros."""
    # Initialize the left-right/right-left data
    X = make_test_data_2(length=length, repeats=1)

    if logits_type == "zeros":
        # Make zero logits
        logits = np.zeros((length, length))
    elif logits_type == "optimal":
        # Make zero logits
        logits = make_left_right_right_left_logits(length)

    # Make HyperTraPS object
    hypertraps = HyperTraPS(logits=logits, verbose=False, n_samples=200)
    actual = hypertraps.log_prob(X)

    # Expected log likelihoods
    if logits_type == "zeros":
        expected = np.array(
            2
            * [
                np.sum(np.log(np.arange(1, i + 1)))
                - np.sum(np.log(np.arange(length, length - i, -1)))
                for i in range(1, length + 1)
            ]
        )[:, np.newaxis]
        assert np.allclose(actual, expected)
    elif logits_type == "optimal":
        expected = np.array([np.log(1 / 2) * (length - 1) * 2])
        assert np.isclose(np.sum(actual), expected)
