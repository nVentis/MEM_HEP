import torch

def kl_div(true:torch.Tensor, test:torch.Tensor, logp:torch.Tensor, logq:torch.Tensor):
    true = true.detach()
    test = test.detach()
    logp = logp.detach()
    
    return (true/test*( logp - logq )).mean()

def exp_div(true:torch.Tensor, test:torch.Tensor, logp:torch.Tensor, logq:torch.Tensor):
    true = true.detach()
    test = test.detach()
    logp = logp.detach()
    
    fac1 = true/test
    fac2 = ( logp - logq )**2
    
    print(fac1.max().item(), fac2.max().item())
    
    return (fac1*fac2).mean()

def chi2_div(true:torch.Tensor, test:torch.Tensor, logp:torch.Tensor, logq:torch.Tensor):
    del logp, logq
    
    return ((true.detach() - test)**2 / test / test.detach()).mean()

def var_loss(true:torch.Tensor, test:torch.Tensor, logp:torch.Tensor, logq:torch.Tensor):
    del logp, logq
    
    return torch.var(true.detach()/test)