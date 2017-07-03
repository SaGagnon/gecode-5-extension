# %load /home/sam/gecode-5.0.0-extension/sql/utils.py
def get_node_in_exs(exs, path_db, sat=True):
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
                index_col=['exec_id','node_id','var_idx','val']
            )
        )

    return df

# REPETITION
features_subset = [
 ('dens', 'maxsd[idx]'),
 ('a_avg_sd', 'aAvgSD[idx]'),
 ('var_dom_size', 'var_dom_size[idx]'),
 ('max_rel_sd', 'maxRelSD[idx]'),
 ('max_rel_ratio', 'maxRelRatio[idx]'),
 ('w_sc_avg', 'wSCAvg[idx]'),
 ('w_anti_sc_avg', 'wAntiSCAvg[idx]'),
 ('w_t_avg', 'wTAvg[idx]'),
 ('w_anti_t_avg', 'wAntiTAvg[idx]'),
 ('w_d_avg', 'wDAvg[idx]')
]
# REPETITION