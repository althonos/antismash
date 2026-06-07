# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running hmmscan.
"""

import io
import logging
from io import StringIO
from typing import List

import pyhmmer

from .base import execute, get_config, SearchIO


def _find_error(output: list[str]) -> str:
    """ Returns the most descriptive line in error output from hmmscan """
    # is there a line that explicitly starts with the error logging?
    for i, line in enumerate(output):
        if line.startswith("Error:"):
            if i + 1 < len(output):
                return f"{line.strip()} {output[i + 1].strip()}"
            return line.strip()
    # if not, take the first non-empty line
    for line in output:
        line = line.strip()
        if line:
            return line
    # in the worst case, return a default
    return "unknown error"


_HMM_CACHE = {}
def _load_hmms(path: str) -> list[pyhmmer.plan7.HMM]:
    if path not in _HMM_CACHE:
        with pyhmmer.plan7.HMMFile(path) as hmm_file:
            _HMM_CACHE[path] = list(hmm_file)
    return _HMM_CACHE[path]


def run_hmmscan(target_hmmfile: str, query_sequence: str, opts: List[str] = None,
                results_file: str = None) -> list[SearchIO._model.query.QueryResult]:
    """ Runs hmmscan on the inputs and return a list of QueryResults

        Arguments:
            target_hmmfile: the path to a HMM file to use in scanning
            query_sequence: a string containing input sequences in fasta format
            opts: a list of extra arguments to pass to hmmscan, or None
            results_file: a path to keep a copy of hmmscan results in, if provided

        Returns:
            a list of QueryResults as parsed from hmmscan output by SearchIO

    """
    if not query_sequence:
        raise ValueError("Cannot run hmmscan on empty sequence")

    logger = logging.getLogger()
    config = get_config()
    cpus = config.cpus

    # Parse query sequences
    query_buffer = io.BytesIO(query_sequence.encode('utf-8'))
    aa = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(
        query_buffer,
        format="fasta",
        digital=True,
        alphabet=aa,
    ) as sequence_file:
        queries = sequence_file.read_block()

    # Pre-load HMMs
    hmms = _load_hmms(target_hmmfile)
    
    # Additional option parsing
    pyhmmer_options = dict(bias_filter=False, cpus=cpus, Z=len(hmms))
    if opts:
        if "--cut_tc" in opts:
            pyhmmer_options["bit_cutoffs"] = "trusted"
            opts.remove("--cut_tc")
        if "--cut_ga" in opts:
            pyhmmer_options["bit_cutoffs"] = "gathering"
            opts.remove("--cut_ga")
        if "-E" in opts:
            i = opts.index("-E")
            pyhmmer_options["E"] = float(opts.pop(i+1))
            opts.pop(i)
        if opts:
            logger.warning(
                "unknown options in run_hmmscan: {}".format(opts)
            )

    # Run hmmscan
    output = io.BytesIO()
    for i, hits in enumerate(pyhmmer.hmmsearch(hmms, queries, **pyhmmer_options)):
        hits.write(output, format="domains", header=i==0)

    # Parse result table
    output.seek(0)
    return list(SearchIO.parse(io.TextIOWrapper(output), "hmmsearch3-domtab"))


def run_hmmscan_help() -> str:
    """ Get the help output of hmmscan """
    # cache results
    help_text = getattr(run_hmmscan_help, 'help_text', '')
    if help_text:
        return help_text

    hmmscan = get_config().executables.hmmscan
    command = [
        hmmscan,
        "-h",
    ]

    help_text = execute(command).stdout
    if not help_text.startswith("# hmmscan"):
        raise RuntimeError(f"unexpected output from hmmscan: {hmmscan}, check path")

    setattr(run_hmmscan_help, 'help_text', help_text)
    return help_text


def run_hmmscan_version() -> str:
    """ Get the version of the hmmscan """
    version_line = run_hmmscan_help().split('\n')[1]
    return version_line.split()[2]
