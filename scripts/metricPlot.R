library(ggplot2)
library(ggpubr)
library(RColorBrewer)


read.tabular <- function(metrics.table.path){

    as.data.frame(read.table(
        file = as.character(metrics.table.path), sep = '\t', header = TRUE))

}

apply.plot.theme <- function(plot, remove.legend=FALSE){

    plot.theme <- plot + theme_pubr() + 
                  theme(text = element_text(size=12, face='bold')) +
                  theme(axis.text.x = element_text(angle = 45, hjust=1))
    if (remove.legend == TRUE){
        plot.theme <- plot.theme + theme(legend.position='none')
    }
    plot.theme


}

plot.colors <- function(metrics.df){
    # https://www.datanovia.com/en/blog/easy-way-to-expand-color-palettes-in-r/
    nb.cols <- nrow(metrics.df)
    colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

}


plot.alignment.lengths <- function(metrics.df, colors){

    
    ggplot(metrics.df, aes(x=read_name, y=length, fill=read_name)) + 
        geom_bar(color='black', size=1, stat='identity') +
        scale_fill_manual(values=colors) +
        theme(legend.position='none') +
        labs(x='Read', y='Length of alignment', title='Length of BLAST alignment') +
        geom_hline(yintercept=280, color='black', size=2)
    
}

plot.EcoRI.SacI.digest.expectation <- function(metrics.df, colors, expected.frag.length=223){


    metrics.df$frag.length.diff <- metrics.df$EcoRI_SacI_digestion_frag_length - expected.frag.length
    metrics.df.primer <- subset(metrics.df, primer=='pFC9_t7_primer_1')
    ggplot(metrics.df, aes(x=read_name, y=frag.length.diff)) + 
        geom_segment(
            aes(x=read_name, xend=read_name, y=0, yend=frag.length.diff), color='black'
            ) +
        geom_point(aes(fill=read_name), colour="black",pch=21, size=5) +
        theme(legend.position='none') +
        labs(
            x='Read', y='Difference in length from expected fragment',
            title='Simulated EcoRI SacI digest')

}

plot.insert.start.distance.expectation <- function(metrics.df, colors){

    metrics.df$distance.expectation <- metrics.df$sstart - metrics.df$expected_insert_distance
    ggplot(metrics.df, aes(x=read_name, y=distance.expectation)) +
        geom_segment(
                aes(x=read_name, xend=read_name, y=0, yend=distance.expectation), color='black'
                ) +
            geom_point(aes(fill=read_name), colour="black",pch=21, size=5) +
            theme(legend.position='none') +
            labs(
                x='Read', y='Difference in expect alignment start position',
                title='Expected alignment start positions')
}


product.plot <- function(df, colors){

    df$ymin <- seq(from = 0.5, to = nrow(df)*2, length.out = nrow(df))
    df$ymax <- df$ymin + 0.75
    print(df$read_length)
    x_max = max(df$read_length)
    template_name = unique(df$qseqid)

    ggplot(df, 
        aes(
            ymin=ymin,
            ymax=ymax,
            xmin=sstart,
            xmax=send,
            fill=read_name)
        ) + xlim(0, x_max) + geom_rect(color='black') + theme_pubr() +
            theme(
                axis.text.y = element_blank(),
                axis.ticks = element_blank()
            ) +
        scale_fill_manual(values=colors) + 
        labs(
            x='Insert reference', 
            y='Sanger sequence alignments', 
            title='Alignment location', 
            fill='Read name')

}

phred.score.boxplot <- function(phred.df, blast.df, colors){


    read.name.df <- blast.df[, c('read_name', 'hash')]
    phred.df.names <- merge(phred.df, read.name.df, by='hash')
    ggplot(phred.df.names, aes(x=read_name, y=phred, fill=read_name)) +
        geom_boxplot(color='black')  + scale_fill_manual(values=colors) +
        labs(title='Mean Phred score', x='', y='Phred score')

}

phred.score.smooth <- function(phred.df, blast.df, colors){


    read.name.df <- blast.df[, c('read_name', 'hash')]
    phred.df.names <- merge(phred.df, read.name.df, by='hash')
    ggplot(phred.df.names, aes(x=position, y=phred, color=read_name)) +
        geom_smooth()  + scale_color_manual(values=colors)

}






save.plot <- function(plot, output.path){


    ggsave(output.path, plot, dpi=600, width=14, height=10, unit='in')


}


main <- function(){


    metrics.filepath <- snakemake@input$metrics_table
    phred.filepath <- snakemake@input$phred_table
    phred.df <- read.tabular(phred.filepath)
    metrics.df <- read.tabular(metrics.filepath)
    if (nrow(metrics.df) > 0){
        colors <- plot.colors(metrics.df)

        # Make all plots
        align.plot <- plot.alignment.lengths(metrics.df, colors)
        digest.plot <- plot.EcoRI.SacI.digest.expectation(metrics.df, colors, 232)
        expected.start.plot <- plot.insert.start.distance.expectation(metrics.df, colors)
        align.plot.vis <- product.plot(metrics.df, colors)
        phred.box <- phred.score.boxplot(phred.df, metrics.df, colors)
        phred.smooth <- phred.score.smooth(phred.df, metrics.df, colors)
        
        # Apply common plot themes 
        align.plot.theme <- apply.plot.theme(align.plot, TRUE)
        digest.plot.theme <- apply.plot.theme(digest.plot, TRUE)
        expected.start.plot.theme <- apply.plot.theme(expected.start.plot, TRUE)
        phred.box.theme <- apply.plot.theme(phred.box, TRUE)
        phred.smooth.theme <-  apply.plot.theme(phred.smooth)
        align.plot.vis.thmeme <- apply.plot.theme(align.plot.vis)
        

        main.plot <- ggarrange(

            ggarrange(align.plot.theme, digest.plot.theme, expected.start.plot.theme, nrow=1, ncol=3),  # top row
            ggarrange(align.plot.vis, ggarrange(phred.box.theme, phred.smooth.theme, nrow=1, ncol=2), nrow=1, ncol=2),  # bottom row
            nrow=2, ncol=1
        )
    
    }else{
        # no alignments and therefore no data to plot
        # is this what we want to do? Maybe have all blast alignments
        # and metrics available in case? For now save empty plot.
        main.plot <- ggplot()
    }
    
    output.pdf <- snakemake@output$pdf
    output.png <- snakemake@output$png

    save.plot(main.plot, output.pdf)
    save.plot(main.plot, output.png)

}

if (sys.nframe() == 0){
    main()
}
