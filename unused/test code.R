# test code for package

# load heartstep data
library(MRTAnalysis)
data(data_mimicHeartSteps)
head(data_mimicHeartSteps)

library(devtools)
load_all()

ID = 'userid'
Ht = c('logstep_30min_lag1','logstep_pre30min','is_at_home_or_work', 'day_in_study')
St = c('logstep_30min_lag1','logstep_pre30min','is_at_home_or_work', 'day_in_study')
At = 'intervention'
outcome = 'logstep_30min'
prob = 'rand_prob'

set.seed(100)
ps = pseudo_outcome_generator_CVlasso(fold = 5, ID, data = data_mimicHeartSteps,
                                      Ht, St, At, prob=prob, outcome, core_num = 5)
my_formula = as.formula(paste("yDR ~ ", paste(St, collapse = " + ")))

library(reticulate)
use_virtualenv("C:/Users/23300/selective-inference/env3")

# Set Python seed
py_run_string("
import random
import numpy as np

random.seed(1)
np.random.seed(1)
")
select = variable_selection_PY_penal_int(ps, ID, my_formula,
                                         virtualenv_path= "C:/Users/23300/selective-inference/env3")
# sometimes the Python can fail which force our R code stops and give error
# the above seed is set for python to reproduce the error

py_run_string("
import random
import numpy as np

random.seed(50)
np.random.seed(50)
")
select = variable_selection_PY_penal_int(ps, ID, my_formula,
                                         virtualenv_path= "C:/Users/23300/selective-inference/env3")


AsyNormbeta_shared = joint_dist_Penal_Int_shared(E = select$E, NE = select$NE,
                                                 pes_outcome = "yDR", data = ps, id = ID, time = "decision_point")
PQR_shared = PQR_Pint_shared(AsyNormbeta_shared, select)

ej = c(0,0,0,1)
AsyNormbeta_ej = joint_dist_Penal_Int_ej(AsyNormbeta_shared, ej)
PQR_ej = PQR_Pint_ej(PQR_shared, AsyNormbeta_ej, AsyNormbeta_shared)
condition = conditional_dist(PQR_shared, PQR_ej, AsyNormbeta_shared, AsyNormbeta_ej, select)

library(devtools)
load_all()

# this seed selection gives error when use Python to find CI
set.seed(100)
py_run_string("
import random
import numpy as np

random.seed(1)
np.random.seed(1)
")
UI_return = DR_WCLS_LASSO(data = data_mimicHeartSteps,
                          fold = 5, ID = ID,
                          time = "decision_point",
                          Ht = Ht, St = St, At = At,
                          prob = prob, outcome = outcome,
                          method_pesu = "CVLASSO",
                          virtualenv_path = "C:/Users/23300/selective-inference/env3",
                          varSelect_program = "Python",
                          standardize_x = T, standardize_y = T)

########################################################
# try other data
# try IHS
df_IHS = read.csv('C:/Users/23300/OneDrive/Desktop/Research/MRTAnalysis/Post-Selective-Inference-Lab/My test code/dfEventsThreeMonthsTimezones.csv')
df_IHS$prob = rep(0.5, length(df_IHS$PARTICIPANTIDENTIFIER))

df_IHS$logStepsPastDay = log(df_IHS$StepsPastDay + 0.5)
df_IHS$logMaxHRPast24Hr = log(df_IHS$maxHRPast24Hours + 0.5)
df_IHS$less27 = (df_IHS$Age < 27)
df_IHS$logHoursSinceLastMood = log(df_IHS$HoursSinceLastMood + 0.5)

##########
# Do log transfer here because those values are too large compared with the outcome
# (1 order of magnitude greater than outcome)
# Otherwise might have singular in: Error in solve.default(t(Q) %*% solve(omega) %*% Q)


ID = 'PARTICIPANTIDENTIFIER'
Ht = c('logStepsPastDay', 'is_weekend', 'logMaxHRPast24Hr', 'logHoursSinceLastMood','less27', "USBorn", "Sex",
       "has_child")
St = c('logStepsPastDay', 'is_weekend', 'logMaxHRPast24Hr',  'logHoursSinceLastMood', 'less27', "USBorn", "Sex",
       "has_child")
At = 'sent'
outcome = 'LogStepsReward'
prob = 'prob'

# remove all missing value in data
df_IHS2 = na.omit(df_IHS[,c(ID, Ht, At, outcome, prob, "time")])

set.seed(100)
py_run_string("
import random
import numpy as np

random.seed(1)
np.random.seed(1)
")
IHS_dataprocessed = DR_WCLS_LASSO(data = df_IHS2,
                          fold = 5, ID = ID,
                          time = "time",
                          Ht = Ht, St = St, At = At,
                          prob = prob, outcome = outcome,
                          method_pesu = "CVLASSO",
                          virtualenv_path = "C:/Users/23300/selective-inference/env3",
                          varSelect_program = "Python")

set.seed(100)
py_run_string("
import random
import numpy as np

random.seed(1)
np.random.seed(1)
")
IHS_dataprocessed_WOstand = DR_WCLS_LASSO(data = df_IHS2,
                                  fold = 5, ID = ID,
                                  time = "time",
                                  Ht = Ht, St = St, At = At,
                                  prob = prob, outcome = outcome,
                                  method_pesu = "CVLASSO",
                                  virtualenv_path = "C:/Users/23300/selective-inference/env3",
                                  varSelect_program = "Python",
                                  standardize_x = F, standardize_y = F)

#######################
# standardization will give differnet confidence interval

set.seed(100)
ps_IHS = pseudo_outcome_generator_CVlasso(fold = 5, ID, data = df_IHS2,
                                      Ht, St, At, prob=prob, outcome, core_num = 5)
my_formula = as.formula(paste("yDR ~ ", paste(St, collapse = " + ")))
select_IHS = variable_selection_PY_penal_int(ps_IHS, ID, my_formula,
                                         virtualenv_path= "C:/Users/23300/selective-inference/env3")

AsyNormbeta_shared_IHS = joint_dist_Penal_Int_shared(E = select_IHS$E, NE = select_IHS$NE,
                                                 pes_outcome = "yDR", data = ps_IHS, id = ID, time = "time",
                                                 moderator_formula = my_formula)
PQR_shared_IHS = PQR_Pint_shared(AsyNormbeta_shared_IHS, select_IHS)

ej_IHS = c(0,1)
AsyNormbeta_ej_IHS = joint_dist_Penal_Int_ej(AsyNormbeta_shared_IHS, ej_IHS)
PQR_ej_IHS = PQR_Pint_ej(PQR_shared_IHS, AsyNormbeta_ej_IHS, AsyNormbeta_shared_IHS)
condition_IHS = conditional_dist(PQR_shared_IHS, PQR_ej_IHS, AsyNormbeta_shared_IHS, AsyNormbeta_ej_IHS,
                                 select_IHS)
####################################################
# check standardization
sim_data = generate_dataset(N = 1000, T = 40, P = 50, sigma_residual = 1.5, sigma_randint = 1.5, main_rand = 3, rho = 0.7,
                                             beta_logit = c(-1, 1.6 * rep(1/50, 50)), model = ~ state1 + state2 + state3 + state4,
                                             beta = matrix(c(-1, 1.7, 1.5, -1.3, -1),ncol = 1),
                                             theta1 = 0.8)
Ht = unlist(lapply(1:50, FUN = function(X) paste0("state",X)))
St = unlist(lapply(1:25, FUN = function(X) paste0("state",X)))
set.seed(100)
py_run_string("
import random
import numpy as np

random.seed(1)
np.random.seed(1)
")

UI_return_sim_notSand = DR_WCLS_LASSO(data = sim_data,
                          fold = 5, ID = "id",
                          time = "decision_point",
                          Ht = Ht, St = St, At = "action",
                          prob = "prob", outcome = "outcome",
                          method_pesu = "CVLASSO",
                          virtualenv_path = "C:/Users/23300/selective-inference/env3",
                          varSelect_program = "Python",
                          standardize_x = F, standardize_y = F)

UI_return_sim_notSand_R = DR_WCLS_LASSO(data = sim_data,
                                      fold = 5, ID = "id",
                                      time = "decision_point",
                                      Ht = Ht, St = St, At = "action",
                                      prob = "prob", outcome = "outcome",
                                      method_pesu = "CVLASSO",
                                      virtualenv_path = "C:/Users/23300/selective-inference/env3",
                                      varSelect_program = "R",
                                      standardize_x = F, standardize_y = F,
                                      beta = matrix(c(-1, 1.7, 1.5, -1.3, -1, rep(0, 21))))

set.seed(100)
py_run_string("
import random
import numpy as np

random.seed(1)
np.random.seed(1)
")
UI_return_sim_Sand = DR_WCLS_LASSO(data = sim_data,
                                      fold = 5, ID = "id",
                                      time = "decision_point",
                                      Ht = Ht, St = St, At = "action",
                                      prob = "prob", outcome = "outcome",
                                      method_pesu = "CVLASSO",
                                      virtualenv_path = "C:/Users/23300/selective-inference/env3",
                                      varSelect_program = "Python",
                                      standardize_x = T, standardize_y = T,
                                   beta = matrix(c(-1, 1.7, 1.5, -1.3, -1, rep(0, 21)),ncol = 1))


# try IHS
ID = 'PARTICIPANTIDENTIFIER'
Ht = c('StepsPastDay', 'is_weekend', 'maxHRPast24Hours', 'HoursSinceLastMood')
St = c('StepsPastDay', 'is_weekend', 'maxHRPast24Hours',  'HoursSinceLastMood')
At = 'sent'
outcome = 'LogStepsReward'
prob = 'prob'


set.seed(200)
py_run_string("
import random
import numpy as np

random.seed(1)
np.random.seed(1)
")
UI_return_IHS = DR_WCLS_LASSO(data = df_IHS,
                              fold = 5, ID = ID,
                              time = "time", Ht = Ht, St = St, At = At,
                              prob = prob, outcome = outcome,
                              virtualenv_path = "C:/Users/23300/selective-inference/env3",
                              #virtualenv_path = "/Users/yuxuanchen/Library/CloudStorage/OneDrive-Personal/Desktop/Research/MRTAnalysis/selective-inference/env3",
                              method_pesu = "CVLASSO", lam = NULL,
                              noise_scale = NULL, splitrat = 0.8,
                              level = 0.9, core_num=3, CI_algorithm = 'lapply',
                              max_iterate = 10^{6}, max_tol = 10^{-3}, varSelect_program = "Python") # default is standardize
# no error returned
set.seed(200)
py_run_string("
import random
import numpy as np

random.seed(1)
np.random.seed(1)
")
notStand_IHS = DR_WCLS_LASSO(data = df_IHS,
                              fold = 5, ID = ID,
                              time = "time", Ht = Ht, St = St, At = At,
                              prob = prob, outcome = outcome,
                              virtualenv_path = "C:/Users/23300/selective-inference/env3",
                              #virtualenv_path = "/Users/yuxuanchen/Library/CloudStorage/OneDrive-Personal/Desktop/Research/MRTAnalysis/selective-inference/env3",
                              method_pesu = "CVLASSO", lam = NULL,
                              noise_scale = NULL, splitrat = 0.8,
                              level = 0.9, core_num=3, CI_algorithm = 'lapply',
                              max_iterate = 10^{6}, max_tol = 10^{-3}, varSelect_program = "Python",
                              standardize_x = F, standardize_y = F) # given error

#######################################################
set.seed(50)
sim_data = generate_dataset_test(N = 100, T = 40, P = 10, sigma_residual = 1.5, sigma_randint = 1.5, main_rand = 3, rho = 0.7,
                            beta_logit = c(-1, 1.6 * rep(1/2000, 5), 1.6 * rep(1/20, 5)), model = ~ state1 + state6,
                            beta = matrix(c(-0.5, 0.01, 1),ncol = 1),
                            theta1 = 0.8)
Ht = unlist(lapply(1:10, FUN = function(X) paste0("state",X)))
St = unlist(lapply(1:10, FUN = function(X) paste0("state",X)))

beta_truth = matrix(c(-0.5, 0.01, rep(0, 4),
                      1, rep(0,4)),ncol = 1)
ps = pseudo_outcome_generator_CVlasso(fold = 5, ID = "id",
                                      data = sim_data, Ht = Ht, St = St, At = "action",
                                      prob="prob", outcome = "outcome", core_num = 5)
hist(ps$prob)
hist(ps$outcome)

sum(is.na(ps$yDR))

library(reticulate)
use_virtualenv("C:/Users/23300/selective-inference/env3")

set.seed(10)
py_run_string("
import random
import numpy as np

random.seed(1)
np.random.seed(1)
")

notSand_test = DR_WCLS_LASSO(data = sim_data,
                                      fold = 5, ID = "id",
                                      time = "decision_point", core_num = 5,
                                      Ht = Ht, St = St, At = "action",
                                      prob = "prob", outcome = "outcome",
                                      method_pesu = "CVLASSO",
                                      virtualenv_path = "C:/Users/23300/selective-inference/env3",
                                      #virtualenv_path = "/Users/yuxuanchen/Library/CloudStorage/OneDrive-Personal/Desktop/Research/MRTAnalysis/selective-inference/env3",
                                      varSelect_program = "Python",
                                      standardize_x = F, standardize_y = F,
                             beta = beta_truth)

mean(notSand_test$post_true <= notSand_test$upperCI & notSand_test$post_true >= notSand_test$lowCI)

set.seed(10)
py_run_string("
import random
import numpy as np

random.seed(1)
np.random.seed(1)
")

Sand_test = DR_WCLS_LASSO(data = sim_data,
                          fold = 5, ID = "id",
                          time = "decision_point", core_num = 5, lam = 0.5,
                          Ht = Ht, St = St, At = "action",
                          prob = "prob", outcome = "outcome",
                          method_pesu = "CVLASSO",
                          virtualenv_path = "C:/Users/23300/selective-inference/env3",
                          #virtualenv_path = "/Users/yuxuanchen/Library/CloudStorage/OneDrive-Personal/Desktop/Research/MRTAnalysis/selective-inference/env3",
                          varSelect_program = "Python",
                          standardize_x = T, standardize_y = T,
                          beta = beta_truth)

mean(Sand_test$post_true <= Sand_test$upperCI & Sand_test$post_true >= Sand_test$lowCI)

# the standardization makes the upper and lower bound very large
# also I have to give very small lam to make sure the algorithm select some variables


# you can also use glment to check the performance of standardization
# use they as benchmark
