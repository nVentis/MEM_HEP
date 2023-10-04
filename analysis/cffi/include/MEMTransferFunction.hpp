class MEMTransferFunction
{
    protected:
        double* params;
        
    public:
        MEMTransferFunction(double pars[]);

        double calc_tf(double reco_val, double true_val);
    };

MEMTransferFunction::MEMTransferFunction(double pars[]) {
    params = pars;
}