# (c) 2014 Brendan Bulik-Sullivan and Hilary Finucane
# Derived from legacy LDSC ldscore/parse.py; GPLv3 (LICENSE).
# Modified 2024 for Python 3; 2026-04 for package extraction and docstrings,
# and 2026-08 for the reduced parser surface. Notice restored 2026-09-10.
# 2026-11: loj selects rows by membership, so a repeated ID counts once;
# .fam and keep-file IDs are read as exact text, without missing values.
# 2026-11-06: the reader class factory became one function returning a
# dataclass; the unused row count ``n`` was dropped.
# See NOTICE for attribution and dated modifications.

"""
plink_files.py

Core functionality:
    Parse PLINK metadata and identifier-list formats used by low-level kernels.

Overview
--------
This module retains only PLINK ``.bim``/``.fam`` and small identifier-list
primitives still required by reference-panel and LD-score computation. Legacy
summary-statistics and LD-score suite parsing belongs to public workflow modules.

Design Notes
------------
- Public regression never imports legacy LD-score/count/annotation fragments.
- The explicit converter, ``ldsc._ldscore.legacy_converter``, owns those compatibility formats.
"""

from dataclasses import dataclass
from functools import partial

import numpy as np
import pandas as pd

from ..errors import LDSCInputError


def get_compression(fh):
    '''Which sort of compression should we use with read_csv?'''
    if fh.endswith('gz'):
        compression = 'gzip'
    elif fh.endswith('bz2'):
        compression = 'bz2'
    else:
        compression = None

    return compression


@dataclass(frozen=True, eq=False)
class IdentifierTable:
    """One identifier-list file, read as legacy LDSC ``parse.py`` read it.

    ``fname`` is the path read, ``df`` holds the selected columns under their
    names and ``IDList`` holds the identifier column as a one-column object
    table, which :meth:`loj` matches against an external list.
    """

    fname: str
    df: pd.DataFrame
    IDList: pd.DataFrame

    def loj(self, externalDf):
        '''Returns indices of those elements of self.IDList that appear in externalDf.

        Rows are selected by membership, in ``self.IDList`` order, so an ID
        listed more than once in ``externalDf`` selects its rows once.
        '''
        return np.flatnonzero(self.IDList.iloc[:, 0].isin(externalDf.iloc[:, 0]).to_numpy())


def read_identifier_table(fname, *, colnames, keepcol, fname_end, usecols, as_text=False) -> IdentifierTable:
    """
    Read one whitespace-delimited, headerless identifier-list file.

    ``usecols`` selects the file's columns, which are named ``colnames``;
    column ``keepcol`` of the result is the identifier column. A non-empty
    ``fname_end`` is the required filename suffix. ``as_text=True`` keeps
    identifiers as their exact text, so ``01003`` and ``1003`` stay distinct
    and words such as ``NA`` or ``nan`` stay identifiers rather than missing
    values.
    """
    if fname_end and not fname.endswith(fname_end):
        raise LDSCInputError(
            f"Legacy LDSC identifier-list reader could not open '{fname}' because the filename does "
            f"not end with '{fname_end}'. Most likely the wrong file type was passed to this loader. "
            f"Pass a file whose name ends with '{fname_end}'."
        )
    text = {'dtype': str, 'na_filter': False} if as_text else {}
    df = pd.read_csv(fname, header=None, usecols=usecols, sep=r'\s+', compression=get_compression(fname), **text)
    df.columns = colnames
    return IdentifierTable(fname, df, df.iloc[:, [keepcol]].astype('object'))


PlinkBIMFile = partial(read_identifier_table, colnames=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], keepcol=1,
                       fname_end='.bim', usecols=[0, 1, 2, 3, 4, 5])
PlinkFAMFile = partial(read_identifier_table, colnames=['IID'], keepcol=0, fname_end='.fam', usecols=[1], as_text=True)
FilterFile = partial(read_identifier_table, colnames=['ID'], keepcol=0, fname_end=None, usecols=[0], as_text=True)
