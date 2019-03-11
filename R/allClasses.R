
setClass('Consensus',
         slots = list(
           IPs = 'tbl_df',
           PIPs = 'list',
           classifications = 'tbl_df',
           consensus = 'tbl_df'
         )
)