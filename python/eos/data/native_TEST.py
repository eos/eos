import unittest

import eos
import pypmc
import numpy as np

class PMCSamplerTests(unittest.TestCase):

    def test_evaluate_mixture_pdf(self):
        "Test the evaluation of a mixture PDF, used in the computation of test statistics."

        component_weights = np.array([0.3, 0.7])

        mean0       = np.array ([ 5.0  , 0.01  ])
        covariance0 = np.array([[ 0.01 , 0.003 ],
                                [ 0.003, 0.0025]])

        mean1       = np.array ([-4.0  , 1.0   ])
        covariance1 = np.array([[ 0.1  , 0.    ],
                                [ 0.   , 0.02  ]])

        component_means = [mean0, mean1]
        component_covariances = [covariance0, covariance1]

        target_mixture = pypmc.density.mixture.create_gaussian_mixture(component_means, component_covariances, component_weights)

        self.assertAlmostEqual(
            eos.data.PMCSampler._evaluate_mixture_pdf(target_mixture, np.array([-3.5, 0.8])),
            0.26256727,
            delta=1e-5
        )

if __name__ == '__main__':
    unittest.main(verbosity=5)
