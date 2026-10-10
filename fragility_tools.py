"""Branch B numerics; input is an explicit internal contract, not a GEM file schema.

Do not use vulnerability means to fabricate the missing exceedance probabilities.
"""
import numpy as np


def state_probabilities(exceedance):
    """Rows: buildings; columns: P(DS>=1), ..., P(DS>=4)."""
    x=np.asarray(exceedance,float)
    if x.ndim!=2 or x.shape[1]!=4 or not np.isfinite(x).all():
        raise ValueError('Four finite DS exceedance probabilities per building required')
    if ((x<0)|(x>1)).any() or (np.diff(x,axis=1)>1e-12).any():
        raise ValueError('Exceedance probabilities must be bounded and nonincreasing across DS')
    return np.column_stack([1-x[:,0],-np.diff(x,axis=1),x[:,-1]])


def score_states(probabilities, observed_ds):
    p=np.asarray(probabilities,float);y=np.asarray(observed_ds)
    if p.shape!=(len(y),5) or not np.isfinite(p).all() or (p<0).any() or not np.allclose(p.sum(axis=1),1):
        raise ValueError('Five valid state probabilities per observation required')
    if not np.isin(y,np.arange(5)).all():raise ValueError('DS labels must be integers 0..4')
    y=y.astype(int);obs=np.eye(5)[y]
    actual=p[np.arange(len(y)),y]
    # Exact zero can arise from floating-point tail truncation of the published
    # lognormal curves.  Apply and disclose a conventional probability floor.
    actual_for_log=np.clip(actual,1e-15,1.0)
    return {'n':len(y),'brier':float(np.mean(np.sum((p-obs)**2,axis=1))),
            'ranked_probability_score':float(np.mean(np.sum((np.cumsum(p,axis=1)[:,:4]-np.cumsum(obs,axis=1)[:,:4])**2,axis=1))),
            'negative_log_score':float(-np.log(actual_for_log).mean())}


def propagate_sequence(first_event_states, conditional_transition):
    """p_after2[j] = sum_i p_after1[i] P(DS_after2=j | DS_after1=i, IM2).

    Transition must come from a damaged-state model. Single-event intact fragilities
    are insufficient. Shared latent capacity must be integrated consistently upstream.
    """
    p=np.asarray(first_event_states,float);t=np.asarray(conditional_transition,float)
    if p.ndim!=2 or p.shape[1]!=5 or t.shape!=(len(p),5,5):raise ValueError('Input shapes invalid')
    if not np.isfinite(p).all() or not np.isfinite(t).all() or (p<0).any() or (t<0).any():raise ValueError('Invalid probabilities')
    if not np.allclose(p.sum(axis=1),1) or not np.allclose(t.sum(axis=2),1):raise ValueError('Probabilities must sum to one')
    if np.any(np.tril(t,k=-1)>1e-12):raise ValueError('No healing between events: downward transitions invalid')
    return np.einsum('ni,nij->nj',p,t)
