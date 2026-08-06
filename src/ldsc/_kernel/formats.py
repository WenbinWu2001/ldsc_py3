"""
formats.py

Core functionality:
    Parse PLINK metadata and identifier-list formats used by low-level kernels.

Overview
--------
This module retains only PLINK ``.bim``/``.fam`` and small identifier-list
primitives still required by reference-panel and LD-score computation. Legacy
sumstats and LD-score suite parsing belongs to public workflow modules.

Design Notes
------------
- Public regression never imports LDSC2 LD-score/count/annotation fragments.
- The explicit converter owns those compatibility formats outside ``_kernel``.
"""

from __future__ import division
import numpy as np
import pandas as pd

from ..errors import LDSCInputError


def read_csv(fh, **kwargs):
    """Read a whitespace-delimited LDSC text file with legacy missing-value rules."""
    return pd.read_csv(fh, sep=r'\s+', na_values='.', **kwargs)


def get_compression(fh):
    '''Which sort of compression should we use with read_csv?'''
    if fh.endswith('gz'):
        compression = 'gzip'
    elif fh.endswith('bz2'):
        compression = 'bz2'
    else:
        compression = None

    return compression


def __ID_List_Factory__(colnames, keepcol, fname_end, header=None, usecols=None):
    """
    Build a small file-reader class for one LDSC identifier-list format.

    The returned class mirrors the historical ``parse.py`` behavior used by the
    LD kernel: it reads one file into ``self.df`` and, when requested, exposes
    a one-column ``IDList`` table that can be left-joined against an external
    list to obtain retained row indices.
    """

    class IDContainer(object):
        """Container that loads one identifier-list table and optional `IDList` view."""

        def __init__(self, fname):
            """Store parser settings and load the requested identifier list."""
            self.__usecols__ = usecols
            self.__colnames__ = colnames
            self.__keepcol__ = keepcol
            self.__fname_end__ = fname_end
            self.__header__ = header
            self.__read__(fname)
            self.n = len(self.df)

        def __read__(self, fname):
            """Read one identifier-list file into ``self.df`` and ``self.IDList``."""
            end = self.__fname_end__
            if end and not fname.endswith(end):
                raise LDSCInputError(
                    f"Legacy LDSC identifier-list reader could not open '{fname}' because the filename does "
                    f"not end with '{end}'. Most likely the wrong file type was passed to this loader. "
                    f"Pass a file whose name ends with '{end}'."
                )

            comp = get_compression(fname)
            self.df = pd.read_csv(fname, header=self.__header__, usecols=self.__usecols__,
                                  sep=r'\s+', compression=comp)

            if self.__colnames__:
                self.df.columns = self.__colnames__

            if self.__keepcol__ is not None:
                self.IDList = self.df.iloc[:, [self.__keepcol__]].astype('object')

        def loj(self, externalDf):
            '''Returns indices of those elements of self.IDList that appear in exernalDf.'''
            r = externalDf.columns[0]
            l = self.IDList.columns[0]
            merge_df = externalDf.iloc[:, [0]]
            merge_df['keep'] = True
            z = pd.merge(self.IDList, merge_df, how='left', left_on=l, right_on=r,
                         sort=False)
            ii = (z['keep'] == True).to_numpy()
            return np.nonzero(ii)[0]

    return IDContainer


PlinkBIMFile = __ID_List_Factory__(['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], 1, '.bim', usecols=[0, 1, 2, 3, 4, 5])
PlinkFAMFile = __ID_List_Factory__(['IID'], 0, '.fam', usecols=[1])
FilterFile = __ID_List_Factory__(['ID'], 0, None, usecols=[0])
AnnotFile = __ID_List_Factory__(None, 2, None, header=0, usecols=None)
ThinAnnotFile = __ID_List_Factory__(None, None, None, header=0, usecols=None)
