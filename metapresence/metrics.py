"""
Metrics calculation for metapresence.
"""


def parse_metrics_criteria(cov, ber, fug1, fug2, nreads, criterion, min_reads, max_for_fug, ber_threshold, fug_threshold, unpaired):
    """
    Determine if a genome is present based on the metric criteria.

    Args:
        cov: Coverage value.
        ber: BER value.
        fug1: FUG1 value.
        fug2: FUG2 value.
        nreads: Number of reads.
        criterion: Criterion for FUG evaluation ('all', 'any', 'mean').
        min_reads: Minimum read threshold.
        max_for_fug: Maximum coverage for using FUG.
        ber_threshold: BER threshold.
        fug_threshold: FUG threshold.
        unpaired: Whether reads are unpaired.

    Returns:
        True if the genome is present, None otherwise.
    """
    if nreads < min_reads:
        return None

    if cov > max_for_fug:
        if ber >= ber_threshold:
            return True
    else:
        if criterion == 'all':
            if ber >= ber_threshold and fug1 >= fug_threshold and (fug2 >= fug_threshold if not unpaired else True):
                return True
        elif criterion == 'any':
            if ber >= ber_threshold and (fug1 >= fug_threshold or (fug2 >= fug_threshold if not unpaired else False)):
                return True
        elif criterion == 'mean':
            if ber >= ber_threshold and ((fug1 + fug2) / 2 >= fug_threshold if not unpaired else fug1 >= fug_threshold):
                return True

    return None