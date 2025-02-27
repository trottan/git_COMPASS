import ROOT
import ROOT.TMath as rm
import numpy as np
import math
from functools import lru_cache
#log sum function (the one we are minimizing)
class logsum(object):
    def __init__(self, cosTheta, phi, Phi, eps, n_Data, cosTheta_sim, phi_sim, Phi_sim, eps_sim, n_sim,batch_size=10):
        self.cosTheta = cosTheta
        self.phi = phi
        self.Phi = Phi
        self.eps = eps
        self.n_Data = n_Data

        self.cosTheta_sim = cosTheta_sim
        self.phi_sim = phi_sim
        self.Phi_sim = Phi_sim
        self.eps_sim = eps_sim
        self.n_sim = n_sim
        self.batch_size = batch_size

    #@lru_cache(maxsize=None)  # Enable caching
    #def cached_call(self, x_tuple):
    #    """Helper function for caching, takes x as a tuple"""
    #    x = np.array(x_tuple)  # Convert tuple back to numpy array
    #    return self._compute_logsum(x)

    #def __call__(self, x):
    #    """Convert x to a tuple and use the cached function"""
    #    return self.cached_call(tuple(x))  # Convert numpy array to tuple

    def compute_batch(self, start_idx, end_idx, dataset='data', x=None):
        N = 0.0
        N_u = 0.0
        N_p = 0.0 
        Pb = 1.0
        if dataset == 'data':
            cT_arr = self.cosTheta[start_idx:end_idx]
            ph_arr = self.phi[start_idx:end_idx]
            Ph_arr = self.Phi[start_idx:end_idx]
            eps_arr = self.eps[start_idx:end_idx]
        else:
            cT_arr = self.cosTheta_sim[start_idx:end_idx]
            ph_arr = self.phi_sim[start_idx:end_idx]
            Ph_arr = self.Phi_sim[start_idx:end_idx]
            eps_arr = self.eps_sim[start_idx:end_idx]
        for i in range(len(cT_arr)):
            cT = cT_arr[i]
            ph = ph_arr[i]
            Ph = Ph_arr[i]
            eps = eps_arr[i]

            sin2Theta = 2 * cT * np.sqrt(1 - cT * cT)
            sinTheta2 =  1 - cT * cT
            N_u = (1/2 * (1 - x[0]) + 1/2 * (3 * x[0] - 1)* cT * cT - rm.Sqrt(2) * x[1] * sin2Theta * rm.Cos(ph) - x[2] * sinTheta2 * rm.Cos(2* ph))
            N_u += (-eps * rm.Cos(2 * Ph) * (x[3] * sinTheta2 + x[4] * cT * cT - rm.Sqrt(2) * x[5] * sin2Theta * rm.Cos(ph) - x[6] * sinTheta2 * rm.Cos(2 * ph)))
            N_u += (-eps * rm.Sin(2 * Ph) * (rm.Sqrt(2) * x[7] * sin2Theta * rm.Sin(ph)+ x[8] * sinTheta2 * rm.Sin(2 * ph))) + (rm.Sqrt(2 * eps * (1 + eps))*rm.Cos(Ph)*(x[9] * sinTheta2 + x[10] * cT * cT - rm.Sqrt(2) * x[11] * sin2Theta * rm.Cos(ph)- x[12] * sinTheta2 * rm.Cos(2 * ph)))
            N_u += (rm.Sqrt(2 * eps * (1 + eps)) * rm.Sin(Ph) * (rm.Sqrt(2) * x[13] * sin2Theta * rm.Sin(ph)+ x[14] * sinTheta2* rm.Sin(2 * ph)))

            N_p = (rm.Sqrt(1 - eps * eps) * (rm.Sqrt(2) * x[15] * sin2Theta * rm.Sin(ph)+ x[16] * sinTheta2 * rm.Sin(2 * ph)))
            N_p += (rm.Sqrt(2 * eps * (1 + eps)) * rm.Cos(Ph) * (rm.Sqrt(2) * x[17] * sin2Theta * rm.Sin(ph)+ x[18] * sinTheta2* rm.Sin(2 * ph)))
            N_p += (rm.Sqrt(2 * eps * (1 + eps)) * rm.Sin(Ph) * (x[19] * sinTheta2 + x[20] * cT * cT - rm.Sqrt(2) * x[21] * sin2Theta * rm.Cos(ph) - x[22] * sinTheta2 * rm.Cos(2 * ph)))
            

            if math.isnan(N_u): print(f"N_u is nan at {i}, {cT},{ph},{Ph},{eps},{dataset}")
            #if math.isnan(N_p): print(f"N_p is nan at {i}")
            #print(N_u + Pb * N_p)
            if (N_u + Pb * N_p) <= 0:
                N += 0
            else:
                if dataset == 'data':
                    N += ROOT.TMath.Log(3 / (8 * rm.Pi() * rm.Pi()) * (N_u + Pb * N_p))
                else:
                    N += 3 / (8 * rm.Pi() * rm.Pi()) * (N_u + Pb * N_p)
            #if math.isnan(N): print(f"N is nan at {i}")
        return N



    def __call__(self, x):

        total_N = 0.0
        total_NF = 0.0

        for start in range(0, self.n_Data, self.batch_size):
            end = min(start + self.batch_size, self.n_Data)
            total_N += self.compute_batch(start, end, dataset='data', x=x)
            #if math.isnan(total_N): print(total_N)
        for start in range(0, self.n_sim, self.batch_size):
            end = min(start + self.batch_size, self.n_sim)
            total_NF += self.compute_batch(start, end, dataset='sim', x=x)

        tot = total_N - self.n_Data* ROOT.TMath.Log(total_NF)
        return -tot
  
ROOT.gInterpreter.Declare("""
    ROOT::Math::Functor foo(const std::function<double(double const *)> &x, int num) {return ROOT::Math::Functor(x,num); }
    //ROOT::Math::Functor foo(GlobalChi2 x, int num) { return ROOT::Math::Functor(x,num); }
    
    """)




def Mini(cosTheta,phi,Phi,eps,n_Data,cosTheta_sim,phi_sim,Phi_sim,eps_sim,n_sim):

    minimum = ROOT.Math.Factory.CreateMinimizer("Minuit2", "")

    minimum.SetMaxFunctionCalls(10000)
    minimum.SetMaxIterations(10000)
    minimum.SetTolerance(1e-5)
    minimum.SetErrorDef(0.5)
    minimum.SetPrintLevel(1)

    lsum = logsum(cosTheta,phi,Phi,eps,n_Data,cosTheta_sim,phi_sim,Phi_sim,eps_sim,n_sim)
    logsumFunctor = ROOT.foo(lsum,23)
    #print(lsum)

    minimum.SetFunction(logsumFunctor)

    #step = np.array([0.01, 0.01, 0.01])
    #variable = np.array([0.11, 0.11, -0.11])
    
    variable_temp = [0.4698,0.2457,-0.2459]
    variable_temp += [0.1769,-0.1662,0.453,0.0362]
    variable_temp += [0.0454,-0.0539,0.0532,0.1456,-0.0376,0.0067,0.0019]
    variable_temp += [0.0027,0.0050,-0.0028,-0.0045,0.0203,-0.0300]
    variable_temp += [-0.0120,-0.0162,0.0163]
    variable = np.array(variable_temp)
    
    #step_temp = [0.01,0.01,0.01]
    #step_temp += [0.01,0.01,0.01,0.01]
    #step_temp += [0.01,0.01,0.01,0.01,0.01,0.01,0.01]
    #step_temp += [0.01,0.01,0.01,0.01,0.01,0.01]
    #step_temp += [0.01,0.01,0.01]
    
    step = np.ones(23) * 0.0001#np.array(step_temp)
    
    
    


    minimum.SetVariable(0,"r0400",variable[0], step[0]);
    minimum.SetVariable(1,"Rer0410",variable[1], step[1]);
    minimum.SetVariable(2,"r041-1",variable[2], step[2]);
    
    minimum.SetVariable(3,"r111",variable[3], step[3]);
    minimum.SetVariable(4,"r100",variable[4], step[4]);
    minimum.SetVariable(5,"Rer110",variable[5], step[5]);
    minimum.SetVariable(6,"r11-1",variable[6], step[6]);
    
    minimum.SetVariable(7,"IMr210",variable[7], step[7]);
    minimum.SetVariable(8,"IMr21-1",variable[8], step[8]);
    
    minimum.SetVariable(9,"r511",variable[9], step[9]);
    minimum.SetVariable(10,"r500",variable[10], step[10]);
    minimum.SetVariable(11,"Rer510",variable[11], step[11]);
    minimum.SetVariable(12,"r51-1",variable[12], step[12]);
    
    minimum.SetVariable(13,"IMr610",variable[13], step[13]);
    minimum.SetVariable(14,"IMr61-1",variable[14], step[14]);
    
    minimum.SetVariable(15,"IMr310",variable[15], step[15]);
    minimum.SetVariable(16,"IMr31-1",variable[16], step[16]);
    
    minimum.SetVariable(17,"IMr710",variable[17], step[17]);
    minimum.SetVariable(18,"IMr71-1",variable[18], step[18]);
    
    minimum.SetVariable(19,"r811",variable[19], step[19]);
    minimum.SetVariable(20,"r800",variable[20], step[20]);
    minimum.SetVariable(21,"r810",variable[21], step[21]);
    minimum.SetVariable(22,"r81-1",variable[22], step[22]);




    #minimum.SetVariableLimits(0,-0.5,0.5);
    #minimum.SetVariableLimits(1,0.01,0.4);
    #minimum.SetVariableLimits(2,-0.4,-0.01);

    
    minimum.Minimize()
    minimum.Hesse()

    xs = minimum.X();
    errs = minimum.Errors();
    for i in range(23):
        print(f"Results {i}: ",xs[i]," +/- ", errs[i])



def GetAngles(fileName):
    fname = fileName
    rdf = ROOT.RDataFrame("tree",fname)

    Phi = rdf.AsNumpy(["Phi"])["Phi"]#.tolist()
    phi = rdf.AsNumpy(["phi"])["phi"]#.tolist()
    cosTheta = rdf.AsNumpy(["cosTheta"])["cosTheta"]#.tolist()
    eps = rdf.AsNumpy(["eps"])["eps"]#.tolist()


    nCount = 10#rdf.Count().GetValue()
    return Phi,phi,cosTheta,nCount,eps


Phi,phi,cosTheta,n_Data,eps = GetAngles("/eos/home-n/ntrotta/DVCS_output/rho_2016/Kin*.root")
Phi_sim,phi_sim,cosTheta_sim,n_sim,eps_sim = GetAngles("/eos/home-n/ntrotta/DVCS_output/rho_2016_MC/Kin*.root")
Mini(cosTheta,phi,Phi,eps,n_Data,cosTheta_sim,phi_sim,Phi_sim,eps_sim,n_sim)
print(len(Phi),len(Phi_sim))

