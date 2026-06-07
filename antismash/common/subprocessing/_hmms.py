import pyhmmer

_HMM_CACHE = {}

def load_hmms(path: str) -> list[pyhmmer.plan7.HMM]:
    if path not in _HMM_CACHE:
        with pyhmmer.plan7.HMMFile(path) as hmm_file:
            if hmm_file.is_pressed():
                hmm_file = hmm_file.optimized_profiles()
            _HMM_CACHE[path] = list(hmm_file)
    return _HMM_CACHE[path]