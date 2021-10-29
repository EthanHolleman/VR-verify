library(plyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)


source('scripts/metricPlot.R')

read.tabular <- function(metrics.table.path){

    as.data.frame(read.table(
        file = as.character(metrics.table.path), sep = '\t', header = TRUE))

}

plot.number.reads.insert <- function(metrics.df, insert.names, colors){
    #print(head(metrics.df))
    count.reads <- count(metrics.df, 'qseqid')
    aligned.inserts <- count.reads$qseqid
    print(count.reads)


    no.alignment <- insert.names[!(insert.names %in% aligned.inserts)]
    # add those that are not included
    df.no.alignment <- data.frame(qseqid=no.alignment, freq=0)
    all.inserts <- rbind(count.reads, df.no.alignment)

    print(all.inserts)

    ggplot(all.inserts, aes(x=qseqid, y=freq, fill=qseqid)) +
        geom_bar(color='black', stat='identity') +
        scale_fill_manual(values=colors) +
        labs(title='Number of successful BLAST alignments by insert', x='', y='')

}

plot.total.reads <- function(samples.df, metrics.df){


    total.reads <- nrow(samples.df)
    align.reads <- nrow(metrics.df)
    align.df <- data.frame(
        read.type=c('total reads', 'aligned'),
        num.reads=c(total.reads, align.reads)
    )
    ggplot(align.df, aes(x=read.type, y=num.reads, fill=read.type)) +
        geom_bar(stat='identity', color='black') + 
        scale_fill_manual(values=c('#9f34eb', '#34eb7d')) +
        labs(x='', y='', title='Total vs aligned reads')


}


plot.confidence <- function(metrics.df){

    count.reads <- count(metrics.df, 'qseqid')
    ggplot(metrics.df, aes(x=qseqid, fill=status)) + 
        geom_bar(color='black') + 
        labs(title='Coutns of assigned confidence read contains insert by insert', x='', y='') +
        scale_fill_brewer(palette='Dark2')

}




main <- function(){


    samples.table.path <- snakemake@params$samples_table
    metrics.table <- snakemake@input$metrics_table
    all.inserts.path <- snakemake@input$insert_list
    output <- as.character(snakemake@output)
    
    # read input tables
    metrics.df <- read.tabular(metrics.table)
    samples.df <- read.tabular(samples.table.path)
    all.inserts <- read.tabular(all.inserts.path)

    insert.colors <- plot.colors(all.inserts)

    # make plots
    total.reads.plot <- plot.total.reads(samples.df, metrics.df)
    aligned.reads <- plot.number.reads.insert(
        metrics.df, all.inserts$insert_name,  insert.colors)
    confidence.plot <- plot.confidence(metrics.df)

    # apply common plot themes
    total.reads.plot.theme <- apply.plot.theme(total.reads.plot, TRUE)
    aligned.reads.theme <- apply.plot.theme(aligned.reads, TRUE)
    confidence.plot.theme <- apply.plot.theme(confidence.plot)

    # arrrange plots into one plot
    main.plot <- ggarrange(
        total.reads.plot.theme, aligned.reads.theme, confidence.plot.theme
    )
    save.plot(main.plot, output)





}

if (sys.nframe() == 0){
    main()
}
