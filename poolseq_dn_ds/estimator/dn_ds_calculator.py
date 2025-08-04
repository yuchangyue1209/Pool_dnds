# estimator/dn_ds_calculator.py
def estimate_dn_ds(variants):
    piN, piS = 0.0, 0.0
    countN, countS = 0, 0
    for v in variants:
        if v['type'] == 'nonsyn':
            piN += v['pi']
            countN += 1
        elif v['type'] == 'syn':
            piS += v['pi']
            countS += 1
    ratio = piN / piS if piS > 0 else float('nan')
    return {"piN": piN, "piS": piS, "dN/dS": ratio, "nonsyn_sites": countN, "syn_sites": countS}
