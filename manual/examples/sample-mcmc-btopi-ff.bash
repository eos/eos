eos-sample-mcmc \
    --global-option model CKMScan \
    --global-option form-factors BCL2008 \
    --scan "CKM::abs(V_ub)"          2e-3   5e-3  --prior flat \
    --scan "B->pi::f_+(0)@BCL2008"   0      1     --prior flat \
    --scan "B->pi::b_+^1@BCL2008"  -20    +20     --prior flat \
    --scan "B->pi::b_+^2@BCL2008"  -20    +20     --prior flat \
    --constraint "B->pi::f_+@IKMvD-2014" \
    --constraint "B^0->pi^+lnu::BR@BaBar-2010B" \
    --constraint "B^0->pi^+lnu::BR@Belle-2010A" \
    --constraint "B^0->pi^+lnu::BR@BaBar-2012D" \
    --constraint "B^0->pi^+lnu::BR@Belle-2013A" \
    --prerun-min     500 \
    --prerun-max    7500 \
    --prerun-update  500 \
    --prerun-only        \
    --output /tmp/mcmc_prerun_btopi+ff.hdf5
#    --chunks          10 \
#    --chunk-size    1000 \
