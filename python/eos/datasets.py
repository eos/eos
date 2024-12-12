# Copyright (c) 2023-2024 Danny van Dyk
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

import eos
import os
import numpy as _np
import yaml as _yaml
import copy as _copy
from .data import MixtureDensity
from .deserializable import Deserializable
from .datasets_file_description import DataSetDescription, PublicLikelihoodDescription

from scipy.special import erf
from scipy.stats import chi2


class DataSets:
    """
    Utility class to download public EOS data sets.
    """

    DOWNLOAD_URL_DOI    = 'https://dx.doi.org/{doi}'
    DOWNLOAD_URL_ZENODO = '{url}/files/eos/data-{id}.zip?download=1'
    DOWNLOAD_URL_GITHUB = 'https://github.com/eos/data/archive/refs/tags/{id}.zip'
    UPDATE_URL          = 'https://github.com/eos/data/raw/{ref}/datasets.yaml'


    def __init__(self, storage_directory=None):
        cache_home = os.getenv('XDG_CACHE_HOME', os.path.join(os.getenv('HOME'), '.cache'))
        if not os.path.exists(cache_home):
            eos.warn(f"The current user's cache home directory {cache_home} does not exist; creating it")
            os.makedirs(cache_home)

        self._data_sets = None

        self.storage_directory = os.path.join(cache_home, 'eos') if storage_directory is None else storage_directory
        if not os.path.exists(self.storage_directory):
            eos.info(f"Creating storage directory for datasets: {self.storage_directory}")
            os.makedirs(self.storage_directory)

        if not os.path.exists(os.path.join(self.storage_directory, 'datasets.yaml')):
            self.update()

        if self._data_sets is None:
            with open(os.path.join(self.storage_directory, 'datasets.yaml')) as f:
                data = _yaml.safe_load(f)
                self._data_sets = { k: DataSetDescription.from_dict(**v) for k, v in data.items() }


    def path(self, id):
        """
        Returns the path to a downloaded EOS data set.

        :param id: The identifier of the EOS data set.
        :type id: str
        :returns: The path to the downloaded and unzipped data set.
        :rtype: str
        """
        result = os.path.join(self.storage_directory, id)

        return result


    def exists(self, id):
        """
        Checks whether an EOS data set has been downloaded.

        :param id: The identifier of the EOS data set.
        :type id: str
        :returns: Whether the data set has been downloaded.
        :rtype: bool
        """
        path = os.path.join(self.storage_directory, id)
        return os.path.exists(path) and os.path.isdir(path) and len(os.listdir(path)) > 0


    def datasets(self):
        """
        Returns a generator to iterate over all the available data sets.
        """
        yield from self._data_sets.items()


    def download(self, id):
        """
        Downloads an EOS data set from Zenodo.org or Github.com

        :param id: The identifier of the data set to download.
        :type id: str
        """
        try:
            import requests
        except ImportError:
            raise RuntimeError("eos.DataSets.download() needs to import the 'requests' module. Please install it.")
        import io, zipfile # part of the Python standard library

        idnoversion = id if 'v' not in id else id[:id.rfind('v')]
        if self._data_sets[id].doi is not None:
            doi             = self._data_sets[id].doi
            url             = self.DOWNLOAD_URL_DOI.format(doi=doi)
            expect_redirect = True
        else:
            url             = self.DOWNLOAD_URL_GITHUB.format(id=id)
            expect_redirect = False

        eos.inprogress(f"Downloading EOS data set '{id}' from '{url}' ...")
        r = requests.get(url)
        if not r.ok:
            raise RuntimeError(f"Could not download EOS data set {id} from {url}; status code: {r.status_code}")

        if expect_redirect:
            if not r.history:
                raise RuntimeError(f"Could not download EOS data set {id} from {url}; expected redirect to Zenodo.org did not happen")

            url = self.DOWNLOAD_URL_ZENODO.format(url=r.url, id=id)
            r = requests.get(url)
            if not r.ok:
                raise RuntimeError(f"Could not download EOS data set {id} from {url}; status code: {r.status_code}")

        directory = os.path.join(self.storage_directory, id)
        eos.inprogress(f"... and extracting to '{directory}' ...")
        os.makedirs(directory, exist_ok=True)
        with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
            for zi in zf.infolist():
                oldname = zi.filename
                newname = os.path.join(*(oldname.split('/'))[1:])
                if newname == '':
                    continue

                zi.filename = newname
                targetdir = os.path.join(self.storage_directory, id)
                zf.extract(zi, path=targetdir)
        eos.completed(f"... done")


    def update(self, ref:str='main'):
        """
        Downloads the registry of data sets from the EOS data repository.

        :param revision: The git reference (either a branch name, a tag, or a GIT revision) of the data repository to download the registry from.
        :type revision: str
        """
        try:
            import requests
        except ImportError:
            raise RuntimeError("eos.DataSets.update() needs to import the 'requests' module. Please install it.")

        update_url = self.UPDATE_URL.format(ref=ref)
        eos.inprogress(f"Downloading data sets information from '{update_url}' ...")
        r = requests.get(update_url)
        if not r.ok:
            raise RuntimeError(f"Could not download data sets information from '{update_url}'; status code: {r.status_code}")
        with open(os.path.join(self.storage_directory, 'datasets.yaml'), 'w') as f:
            f.write(r.text)
        eos.completed(f"... done")

        data = _yaml.safe_load(r.text)
        self._data_sets = { k: DataSetDescription.from_dict(**v) for k, v in data.items() }

    def _chi2interp(self, mllh, test_statistics, dof):
        """
        Internal function that evaluate the chi^2 for dof degrees of freedom for a value of negative log-likelihood

        :param mllh: Value of (-2 * log L) for which the chi2 for dof degrees of freedom is returned
        :type mllh: double
        :param test_statistics: Test statistics (i.e. map between pulls in standard deviation and negative log likelihood values) used for the evaluation.
        :type test_statistics: dict
        :param dof: Number of degrees of freedom of the test statistics.
        :type dof: int
        """

        chi2_values = chi2.isf(1 - erf(test_statistics["sigma"] / _np.sqrt(2)), df = dof)
        mllh_values = test_statistics["densities"]

        if mllh < _np.min(mllh_values):
            eos.warn("Found unexpected value for the likelihood when evaluating the test statistics")

        if mllh > _np.max(mllh_values):
            eos.warn(f"Likelihood value exceed the test statictics support. Assigning the maximal chi2 = {max(chi2_values):.2f}")

        return _np.interp(mllh, mllh_values, chi2_values)


    def likelihood(self, id:str, likelihood:str):
        """
        Creates a likelihood function from a public EOS data set.

        :param id: The identifier of the EOS data set.
        :type id: str
        :param likelihood: The name of the likelihood within the data set.
        :type likelihood: str
        """
        if likelihood not in self._data_sets[id].likelihoods:
            raise ValueError(f'Likelihood {likelihood} not found in data set {id}')

        likelihood = self._data_sets[id].likelihoods[likelihood]
        if likelihood.filetype == 'MixtureDensity':
            f = eos.MixtureDensity(os.path.join(self.storage_directory, id, likelihood.filename))
            varied_parameters = f.varied_parameters
            log_pdf = f.density()
            test_statistics = f.test_statistics
            if test_statistics:
                return (
                    varied_parameters,
                    lambda x: -log_pdf.evaluate(x),
                    lambda x: self._chi2interp(x, test_statistics, len(varied_parameters))
                )
            else:
                return (
                    varied_parameters,
                    lambda x: -log_pdf.evaluate(x),
                    None
                )
        elif likelihood.filetype == 'NabuLikelihood':
            import nabu
            f = eos.NabuLikelihood(os.path.join(self.storage_directory, id, likelihood.filename))
            varied_parameters = f.varied_parameters
            log_pdf = f.likelihood.log_prob
            chi2 = f.likelihood.chi2
            return (
                varied_parameters,
                lambda x: -log_pdf(x),
                chi2
            )
        else:
            raise ValueError(f'Unsupported likelihood file type: {likelihood.filetype}')


    def _repr_html_(self):
        result = '<table>\n'
        result += '<tr><th>id</th><th>title</th></tr>\n'
        for id, dataset in self._data_sets.items():
            result += f'<tr><td>{id}</td><td>{dataset.title}</td></tr>\n'
        result += '</table>'

        return result
