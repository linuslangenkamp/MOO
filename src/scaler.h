
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
    std::unique_ptr<double> evalMR_a_scaling; // shift scaling vector
    std::unique_ptr<double> evalMR_D_scaling; // diagonal scaling matrix
};

void NLP::scaleJacobianData() {
    // in place

};

void NLP::scaleHessianData() {
    // in place

};
*/