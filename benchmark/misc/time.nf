Channel
    .fromPath(params.trace)
    .splitCsv(header:true, sep:'\t')
    .map{ row -> row.realtime != '-' ? nextflow.util.Duration.of(row.realtime).millis * (row.cpus as int) : 0 }
    .sum()
    .println { nextflow.util.Duration.of(it).minutes }
