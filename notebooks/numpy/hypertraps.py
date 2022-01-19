"""Class for HyperTraPS inference."""
import numpy as np
from scipy.special import softmax

# Currently not implemented to facilitate numba
# from numba import jit


class HyperTraPS:
    def __init__(self, logits=None, length=None, **kwargs):
        if logits is None and length is None:
            raise ("Logits or length must be specified.")

        self.n_samples = kwargs.get("n_samples", 20)
        self.verbose = kwargs.get("verbose", False)
        self.params_range = kwargs.get("params_range", (-10, 10))

        if logits is not None:
            self.length = logits.shape[-1]
            self.logits = logits

    def _split_data(self, X):
        # Add samples axis to start and end
        start = np.repeat(X[:, 0][np.newaxis, ...], self.n_samples, axis=0)
        end = np.repeat(X[:, 1][np.newaxis, ...], self.n_samples, axis=0)

        return start, end

    def _move_probs(self, off_diagonal, on_diagonal, current, mask_):
        prob_move_per_site = (
            current @ off_diagonal
            + np.diagonal(on_diagonal, axis1=-2, axis2=-1)
        )[0, ..., np.newaxis] * mask_ + np.finfo(
            dtype=self.logits.dtype
        ).minexp * (
            1 - mask_
        )
        return prob_move_per_site

    # @jit(nopython=True)
    def _alphas_loop(self, start, end, on_diagonal, off_diagonal):
        # Alphas matrix
        alphas = []
        current = start.copy()
        for i in range(self.length):
            # Tile the mask across data
            mask = (end != current)[..., np.newaxis]
            mask_all = (np.ones_like(start) != current)[..., np.newaxis]

            # Get prob of moves
            p_all = self._move_probs(
                off_diagonal, on_diagonal, current, mask_all
            )

            # Get weights for compatible moves
            alphas.append((softmax(p_all, axis=2) * mask).sum(axis=2))

            # Get compatible move probs
            p = self._move_probs(off_diagonal, on_diagonal, current, mask)
            p = softmax(p, axis=2)

            # Get uniform random numbers to broadcast over the probs vector
            randoms = np.random.uniform(size=(*p.shape[:2], 1, 1))

            # Calculate where the move occurs
            moves = (randoms < p.cumsum(axis=2)) * 1
            moves = np.concatenate(
                [moves[:, :, :1], np.diff(moves, axis=2)], axis=2
            )[..., 0]
            moves *= np.clip((end != current).sum(axis=2), a_min=0, a_max=1)[
                ..., np.newaxis
            ]

            # Update current positions
            current += moves

        # Convert the alphas into numpy array
        alphas = np.transpose(np.array(alphas), (1, 2, 0, -1))
        return alphas

    def _calculate_alphas(self, X):
        # Get split of data
        start, end = self._split_data(X)

        # Set-up the parameter matrices
        eye = np.eye(*self.logits.shape)
        off_diagonal = self.logits[np.newaxis, np.newaxis, ...] * (1 - eye)
        on_diagonal = self.logits[np.newaxis, np.newaxis, ...] * eye

        # Alphas matrix
        return self._alphas_loop(start, end, on_diagonal, off_diagonal)

    def _alphas_to_log_probs(self, X, alphas, verbose=False):
        # Average over the n_samples dim
        alphas = np.mean(alphas, axis=0)

        # Mask out non-moves
        start_count = np.sum(X[np.newaxis, :, 0], axis=-1)
        end_count = np.sum(X[np.newaxis, :, 1], axis=-1)

        # Get mask for where moves occurred
        alphas_mask = np.ones_like(alphas)
        alphas_mask = np.multiply(
            (start_count.T <= np.arange(self.length)[np.newaxis, ...]),
            (end_count.T > np.arange(self.length)[np.newaxis, ...]),
        )[..., np.newaxis]

        if self.verbose:
            print(alphas)

        # Get log likelihoods at move positions
        alphas[alphas_mask] = np.log(alphas[alphas_mask])
        alphas[~alphas_mask] = np.nan

        # Get sum over moves
        log_probs = np.nansum(alphas, axis=1)

        return log_probs

    def log_prob(self, X, seed=None):
        if seed is not None:
            # Seed
            np.random.seed(seed)

        # Calculate alphas
        alphas = self._calculate_alphas(X)

        # Calculate and return log_probs
        log_probs = self._alphas_to_log_probs(X, alphas)

        return log_probs

    def prob(self, X):
        return np.exp(self.log_prob(X))

    def sample(self, X=None, seed=None):
        # TODO: implement batch samples both with and without conditioning on
        # start and end states optionally passed in X
        pass
