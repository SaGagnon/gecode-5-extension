
#%load /home/sam/gecode-5.0.0-extension/analysis/utils.py
#%%writefile /home/sam/gecode-5.0.0-extension/analysis/utils.py
from sklearn.linear_model import LogisticRegressionCV

import math
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from io import StringIO
%matplotlib inline

def get_node_in_exs(exs, path_db):
    req_sql = """
        SELECT d.*,
          CASE WHEN r.exec_id IS NOT NULL THEN 1 ELSE 0 END as in_sol
        FROM densities AS d
        LEFT JOIN results AS r
          ON d.exec_id=r.exec_id
          AND d.var_id=r.var_id
          AND d.val=r.val
        WHERE d.exec_id = $2;
    """

    df = pd.DataFrame()
    for ex in exs:
        req_sql_ex = req_sql.replace('$2', str(ex))
        output = !sqlite3 -header -csv {path_db} "{req_sql_ex}"
        if len(output) == 0: continue
        df = df.append(
            pd.read_csv(
                StringIO(output.n),
                index_col=['exec_id','node_id','var_id','val']
            )
        )
    return df

def get_node_in_exs_old_db(exs, path_db, sat=True):
    req_sql = """
        SELECT d.*,
          n.sat,
          CASE WHEN r.exec_id IS NOT NULL THEN 1 ELSE 0 END as in_sol
        FROM densities AS d
        JOIN nodes AS n
          ON d.exec_id=n.exec_id
          AND d.node_id=n.node_id
          $1
        LEFT JOIN results AS r
          ON d.exec_id=r.exec_id
          AND d.var_idx=r.var_idx
          AND d.val=r.val
          AND r.res_id=0 -- TEMPORAIRE
        WHERE d.exec_id = $2
           AND EXISTS (
             SELECT exec_id
             FROM results as rr
             WHERE rr.exec_id = $2
        );
    """
    
    if sat:
        req_sql = req_sql.replace('$1',"AND n.sat=1")
        
    df = pd.DataFrame()
    for ex in exs:
        req_sql_ex = req_sql.replace('$2', str(ex))
        output = !sqlite3 -header -csv {path_db} "{req_sql_ex}"
        if len(output) == 0:continue
        df = df.append(
            pd.read_csv(
                StringIO(output.n),
                index_col=['exec_id','node_id','prop_id','var_idx','val']
            )
        )
        
    return df

features_subset = [
    "max_sd",
    "a_avg_sd"
]

def plot_features_sln_sep(features):
    width = 3
    height = math.ceil(len(features)/width)
    plt.figure(figsize=(16,4*height))

    for i, feature in enumerate(features):
        plt.subplot(width,height, i+1)
        plt.title(feature)
        sns.kdeplot(df[df.in_sol == False][feature], color='r')
        sns.kdeplot(df[df.in_sol == True][feature], color='g')
        plt.gca().legend_.remove()
        plt.ylim(0,10)
        
def get_X_y(df):
    return df.iloc[:,:-1], df.iloc[:,-1]

# TEMP
def print_coefs(clf, features):
    print('double _x = 0;')
    for i, coef in enumerate(clf.coef_[0]):
         print("_x += %.4f * %s;" % (coef, features[i][1]))
    print('double intercept = %.4f;' % (clf.intercept_))
    print('_x += intercept;')

# TEMP