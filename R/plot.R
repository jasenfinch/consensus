#' Sankey plot
#' @description
#' Plot an overview of consensus structural classifcations using a sankey plot.
#' @param x The consensus structural classifications to plot. This should either a tibble of an object of S4 class `Construction` as returned by `construction()`.
#' @param exclude_kingdoms a vector of kingdoms names to exclude from the plot
#' @examples
#' plotSankey(construction::structural_classifications)
#' @importFrom ggsankey make_long geom_sankey geom_sankey_text theme_sankey
#' @importFrom ggplot2 ggplot aes theme element_text guides labs
#' @importFrom purrr flatten_chr
#' @importFrom dplyr across
#' @export

setGeneric('plotSankey',function(x,
                                 exclude_kingdoms = c('Unclassified','No database hits','Unknown'))
  standardGeneric('plotSankey')
)

#' @rdname plotSankey

setMethod('plotSankey',
          signature = 'tbl_df',
          function(x,exclude_kingdoms = c('Unclassified','No database hits','Unknown')){
            
            x <- x %>% 
              filter(
                !kingdom %in% exclude_kingdoms 
              ) 
            
            df <- x %>% 
              make_long(
                any_of(c('kingdom',
                       'superclass',
                       'class',
                       'subclass')),
                contains('level')) %>% 
              filter(
                !(is.na(node) & is.na(next_node))
              )
            
            node_orders <- x %>% 
              select(
                any_of(c('kingdom',
                       'superclass',
                       'class',
                       'subclass')),
                contains('level')) %>% 
              distinct() %>% 
              arrange(
                across(
                  everything()
                )
              ) %>% 
              as.list() %>% 
              map(
                ~.x %>% 
                  na.omit() %>% 
                  unique()
              ) %>% 
              flatten_chr()
            
            
            df %>% 
              mutate(
                node = factor(
                  node,
                  levels = node_orders
                ),
                next_node = factor(
                  next_node,
                  levels = node_orders
                )
              ) %>% 
              ggplot(
                aes(x = x, 
                    next_x = next_x, 
                    node = node, 
                    next_node = next_node,
                )
              ) +
              geom_sankey(
                flow.fill = '#B5E6FF',
                node.fill = 'lightgrey'
              ) +
              geom_sankey_text(
                aes(
                  label = node
                ),
                size = 3) +
              theme_sankey() +
              theme(
                axis.title.x = element_text(
                  face = 'bold'
                ),
                axis.text.x = element_text(
                  colour = 'black'
                )
              ) +
              guides(
                fill = 'none'
              ) +
              labs(
                x = 'CHEMONT taxonomic level'
              )
            
          })

#' @rdname plotSankey

setMethod('plotSankey',
          signature = 'Construction',
          function(x,exclude_kingdoms = c('Unclassified','No database hits','Unknown')){
            
            x %>% 
              classifications() %>% 
              plotSankey(exclude_kingdoms = exclude_kingdoms)
          })