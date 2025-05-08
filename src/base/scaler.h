
// TODO: add this later, currently not needed for initial workflow
// should be able to perform at least nominal scaling globally

/*
void NLP::scaleCurrentIterate() {
    // x_scaled = x_D_scaling^(-1) * (x_set_unscaled + (-1) * x_a_scaling)
    Linalg::dsaxpy(numberVars, x_set_unscaled.get(), x_a_scaling.get(), x_D_scaling.get(), -1, true, x_scaled.get());
};

void NLP::unscaleCurrentIterate() { //D * x + beta * y or D^(-1) * x + beta * y
    // x_set_unscaled = x_D_scaling * x_scaled + 1 * x_a_scaling
    Linalg::dgmv(numberVars, x_scaled.get(), x_a_scaling.get(), x_D_scaling.get(), 1, false, x_set_unscaled.get());
};

void NLP::scaleEvalData() {
    // in place
    std::unique_ptr<f64> evalMR_a_scaling; // shift scaling vector
    std::unique_ptr<f64> evalMR_D_scaling; // diagonal scaling matrix
};

void NLP::scaleJacobianData() {
    // in place

};

void NLP::scaleHessianData() {
    // in place

};
*/

    /* TODO: add external scaler class which can perform, no, nominal, adaptive scaling
    // TODO: use these later, fill one time and then scale at the end of calculations
    // these have the same sizes as the curr_'s, just divide element wise
    //FixedVector<f64> curr_x_unscaled;
    f64 curr_obj_nominal = 1;
    FixedVector<f64> curr_grad_nominal;
    FixedVector<f64> curr_g_nominal;
    FixedVector<f64> curr_jac_values_nominal;
    FixedVectorf64> curr_hes_values_nominal;
    */