# test code for package

# load heartstep data
library(MRTAnalysis)
data(data_mimicHeartSteps)
head(data_mimicHeartSteps)

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
set.seed(100)
py_run_string("
import random
import numpy as np

random.seed(50)
np.random.seed(50)
")
UI_return = DR_WCLS_LASSO(data = data_mimicHeartSteps,
                          fold = 5, ID = ID,
                          time = "decision_point",
                          Ht = Ht, St = St, At = At,
                          prob = prob, outcome = outcome,
                          method_pesu = "CVLASSO",
                          virtualenv_path = "C:/Users/23300/selective-inference/env3",
                          varSelect_program = "Python")

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

random.seed(50)
np.random.seed(50)
")
UI_return_IHS = DR_WCLS_LASSO(data = df_IHS2,
                          fold = 5, ID = ID,
                          time = "time",
                          Ht = Ht, St = St, At = At,
                          prob = prob, outcome = outcome,
                          method_pesu = "CVLASSO",
                          virtualenv_path = "C:/Users/23300/selective-inference/env3",
                          varSelect_program = "Python")

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
