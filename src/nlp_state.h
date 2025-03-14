#ifndef OPT_STATE_H
#define OPT_STATE_H

// simple state to check which actions are / have to be performed for an iteration
struct NLP_State {
    bool x_unscaled = false;
    bool eval_f     = false;
    bool eval_g     = false;
    bool grad_f     = false;
    bool jac_g      = false;
    bool hes_lag    = false;

    void check_reset(bool new_x) {
        if (new_x) {
            x_unscaled = false;
            eval_f     = false;
            eval_g     = false;
            grad_f     = false;
            jac_g      = false;
            hes_lag    = false; 
        }
    };
    
    inline bool eval_performed() {
        return (eval_f && eval_g);
    };

    inline bool jac_performed() {
        return (grad_f && jac_g);
    };

    inline bool hes_performed() {
        return hes_lag;
    }
};

#endif  // OPT_STATE_H