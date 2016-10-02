eos-propagate-uncertainty \
    --kinematics s_min  0.0 \
    --kinematics s_max 12.0 \
    --observable "B->pilnu::BR;l=e,form-factors=BCL2008" \
    --pmc-input /tmp/pmc_monolithic_btopi+ff.hdf5 0 10000 \
    --pmc-sample-directory "/data/final/" \
    --output /tmp/unc_btopilnu.hdf5
