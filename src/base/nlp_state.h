#ifndef OPT_STATE_H
#define OPT_STATE_H

// simple state to check which actions are / have to be performed for an iteration
struct NLP_State {
    bool x_set_unscaled = false;
    bool lambda_set     = false;
    bool eval_f         = false;
    bool eval_g         = false;
    bool grad_f         = false;
    bool jac_g          = false;
    bool hes_lag        = false;

    void check_reset_x(bool new_x) {
        if (new_x) {
            x_set_unscaled = false;
            lambda_set     = false;
            eval_f         = false;
            eval_g         = false;
            grad_f         = false;
            jac_g          = false;
            hes_lag        = false; 
        }
    };

    void check_reset_lambda(bool new_lambda) {
        if (new_lambda) {
            lambda_set     = false;
            hes_lag        = false; 
        }
    };
};

#endif  // OPT_STATE_H