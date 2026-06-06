# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running hmmsearch.
"""

import logging
import io
import os
from typing import List

import pyhmmer
from helperlibs.wrappers.io import TemporaryDirectory

from .base import execute, get_config, SearchIO

logger = logging.getLogger()


def run_hmmsearch(query_hmmfile: str, target_sequence: str, use_tempfile: bool = False
                  ) -> List[SearchIO._model.query.QueryResult]:  # pylint: disable=protected-access
    """ Run hmmsearch on a HMM file and a fasta input

        Arguments:
            query_hmmfile: the path to the HMM file
            target_sequence: the fasta input to search as a string
            use_tempfile: if True, a tempfile will be written for the fasta input
                          instead of piping

        Returns:
            a list of hmmsearch results as parsed by SearchIO
    """
    config = get_config()
    cpus = config.cpus

    # Allow for disabling multithreading for HMMer3 calls in the command line
    if config.get('hmmer3') and 'multithreading' in config.hmmer3 and \
            not config.hmmer3.multithreading:
        cpus = 1

    # Parse target sequences
    target_buffer = io.BytesIO(target_sequence.encode('utf-8'))
    aa = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(
        target_buffer,
        format="fasta",
        digital=True,
        alphabet=aa,
    ) as sequence_file:
        targets = sequence_file.read_block()

    # Run hmmsearch
    output = io.BytesIO()
    with pyhmmer.plan7.HMMFile(query_hmmfile) as hmms:
        for i, hits in enumerate(pyhmmer.hmmsearch(hmms, targets, cpus=cpus)):
            hits.write(output, format="domains", header=i==0)

    # Return results
    output.seek(0)
    return list(SearchIO.parse(io.TextIOWrapper(output), "hmmsearch3-domtab"))


def run_hmmsearch_version() -> str:
    """ Get the version of the hmmsearch """

    hmmsearch = get_config().executables.hmmsearch
    command = [
        hmmsearch,
        "-h",
    ]

    help_text = execute(command).stdout
    if not help_text.startswith("# hmmsearch"):
        msg = "unexpected output from hmmsearch: %s, check path"
        raise RuntimeError(msg % hmmsearch)

    version_line = help_text.split('\n')[1]
    return version_line.split()[2]
