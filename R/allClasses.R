
setClass('Consensus',
         slots = list(
           hits = 'tbl_df',
           PIPs = 'tbl_df',
           classifications = 'tbl_df',
           consensus = 'tbl_df'
         )
)

setClass('Consensuses',
         slots = list(
           consensuses = 'list',
           results = 'tbl_df'
         )
)