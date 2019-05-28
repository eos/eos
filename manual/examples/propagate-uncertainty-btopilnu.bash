eos-propagate-uncertainty \
    --kinematics q2  0.00 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  0.25 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  0.50 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  0.75 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  1.00 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  1.50 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  2.00 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  2.50 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  3.00 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  3.50 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  4.00 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  6.00 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  7.00 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  8.00 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2  9.00 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2 10.00 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2 11.00 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --kinematics q2 12.00 --observable "B->pilnu::dBR/dq2;l=e,form-factors=BCL2008,model=CKMScan" \
    --pmc-input /tmp/pmc_monolithic_btopi+ff.hdf5 0 10000 \
    --pmc-sample-directory "/data/final/" \
    --output /tmp/unc_btopilnu.hdf5
