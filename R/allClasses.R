
setClass('Consensus',
         slots = list(
           IPs = 'tbl_df',
           PIPs = 'tbl_df',
           classifications = 'tbl_df',
           consensus = 'list'
         )
)