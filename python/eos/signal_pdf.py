# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2021 Danny van Dyk
#
# This file is part of the EOS project. EOS is free software;
# you can redistribute it and/or modify it under the terms of the GNU General
# Public License version 2, as published by the Free Software Foundation.
#
# EOS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA  02111-1307  USA

from _eos import _SignalPDF, _SignalPDFs
import eos
import numpy as np

class SignalPDF(_SignalPDF):
    def log_pdf(self, x):
        """
        Adapter for use with external optimization software (e.g. pypmc) to aid when sampling from the log(PDF).

        :param x: Phase space point, with the elements corresponding to the sorted list of variables in lexicographical order; inspect `self.variables`.
        :type x: iterable of float
        """
        for kv, v in zip(self.variables, x):
            kv.set(v)

        return self.evaluate()

    def sample_mcmc(self, N, stride, pre_N, preruns, cov_scale=0.1, start_point=None, rng=np.random.mtrand):
        """
        Return samples of the kinematic variables and the log(PDF).

        Obtains random samples of the log(PDF) using an adaptive Markov Chain Monte Carlo with PyPMC.
        A prerun with adaptations is carried out first and its samples are discarded.

        :param N: Number of samples that shall be returned
        :param stride: Stride, i.e., the number by which the actual amount of samples shall be thinned to return N samples.
        :param pre_N: Number of samples in each prerun.
        :param preruns: Number of preruns.
        :param cov_scale: Scale factor for the initial guess of the covariance matrix.
        :param start_point: Optional starting point for the chain
        :type start_point: list-like, optional
        :param rng: Optional random number generator (must be compatible with the requirements of pypmc.sampler.markov_chain.MarkovChain)

        :return: A tuple of the kinematic variables as array of size N and the log(PDF) as array of size N.

        .. note::
           This method requires the PyPMC python module, which can be installed from PyPI.
        """
        import pypmc
        try:
            from tqdm.auto import tqdm
            progressbar = tqdm
        except ImportError:
            progressbar = lambda x, **kw: x

        ind_lower = np.array([bound[0].evaluate() for bound in self.bounds])
        ind_upper = np.array([bound[1].evaluate() for bound in self.bounds])
        ind = pypmc.tools.indicator.hyperrectangle(ind_lower, ind_upper)

        log_target = pypmc.tools.indicator.merge_function_with_indicator(self.log_pdf, ind, -np.inf)

        # create initial covariance
        sigma = np.diag([np.square(bound[1].evaluate() - bound[0].evaluate()) / 12 * cov_scale for bound in self.bounds])
        log_proposal = pypmc.density.gauss.LocalGauss(sigma)

        # create start point, if not provided
        if start_point is None:
            u = np.array([rng.uniform(0.0, 1.0) for j in range(0, len(ind_lower))])
            ubar = 1.0 - u
            start_point = ubar * ind_upper + u * ind_lower

        # create MC sampler
        sampler = pypmc.sampler.markov_chain.AdaptiveMarkovChain(log_target, log_proposal, start_point, save_target_values=True, rng=rng)

        # pre run to adapt markov chains
        for i in progressbar(range(0, preruns), desc="Pre-runs", leave=False):
            eos.info('Prerun {} out of {}'.format(i, preruns))
            accept_count = sampler.run(pre_N)
            accept_rate  = accept_count / pre_N * 100
            eos.info('Prerun {}: acceptance rate is {:3.0f}%'.format(i, accept_rate))
            sampler.adapt()
        sampler.clear()

        # obtain final samples
        eos.info('Main run: started ...')
        sample_total  = N * stride
        sample_chunk  = sample_total // 100
        sample_chunks = [sample_chunk for i in range(0, 99)]
        sample_chunks.append(sample_total - 99 * sample_chunk)
        for current_chunk in progressbar(sample_chunks, desc="Main run", leave=False):
            accept_count = accept_count + sampler.run(current_chunk)
        accept_rate  = accept_count / (N * stride) * 100
        eos.info('Main run: acceptance rate is {:3.0f}%'.format(accept_rate))

        parameter_samples = sampler.samples[:][::stride]
        weights = sampler.target_values[:][::stride, 0]

        return(parameter_samples, weights)

    @staticmethod
    def make(name, parameters, kinematics, options):
        """
        Makes a new :class:`SignalPDF` object.

        :param name: The name of the probability density function (PDF). See `the complete list of PDFs <../signal-pdfs.html>`_.
        :type name: eos.QualifiedName
        :param parameters: The set of parameters to which this PDF is bound.
        :type parameters: eos.Parameters
        :param kinematics: The set of kinematic variables to which this PDF is bound.
        :type kinematics: eos.Kinematics
        :param options: The set of options relevant to this PDF.
        :type options: eos.Options

        :return: The new PDF object.
        :rtype: eos.SignalPDF
        """
        pdf = _SignalPDF.make(name, parameters, kinematics, options)
        pdf.__class__ = SignalPDF
        pdf.variables = list(map(
            lambda n: kinematics[n],
            filter(lambda n: not(n.endswith('_min') or n.endswith('_max')), [kv.name() for kv in kinematics])
        ))
        pdf.bounds = [(kinematics[v.name() + '_min'], kinematics[v.name() + '_max']) for v in pdf.variables]

        return pdf

class SignalPDFs(_SignalPDFs):
    """
    Represents the complete list of probability density functions (PDFs) known to EOS.

    Objects of this class are visualized as tables in Jupyter notebooks for easy
    overview. Filters can be applied through keyword arguments to the initialization.

    :param prefix: Only show observables whose qualified names contain the provided ``prefix`` in their prefix part.
    :type prefix: str
    :param name: Only show observables whose qualified names contain the provided ``name`` in their name part.
    :type name: str
    :param suffix: Only show observables whose qualified names contain the provided ``suffix`` in their suffix part.
    :type suffix: str

    See also `the complete list of signal PDFs <../signal-pdfs.html>`_ in this documentation.
    """
    def __init__(self, prefix=None, name=None, suffix=None, showall=False):
        super().__init__()
        self.prefix=prefix
        self.name=name
        self.suffix=suffix
        self.showall=showall

    def filter_entry(self, qn):
        if self.prefix and not self.prefix in str(qn.prefix_part()):
            return False

        if self.name and not self.name in str(qn.name_part()):
            return False

        if self.suffix and not self.suffix in str(qn.suffix_part()):
            return False

        return True

    def _repr_html_(self):
        result = '<table>'
        for qn, entry in self:
            if not self.filter_entry(qn):
                continue
            result += r'      <tr><th><tt style="color:grey">{qn}</tt></th><td style="text-align:left">{desc}</td></tr>'.format(qn=qn,desc=entry.description())
        result += '</table>'

        return result