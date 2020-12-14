library(EILA)

## Two ancestries
data(ceuyri)
res.eila <- eila(admixed  = ceuyri$admixed,
                 position = ceuyri$position,
                 anc1     = ceuyri$anc1,
                 anc2     = ceuyri$anc2)
cat("Overall accuracy:", mean((res.eila$local.ancestry ==
                                 ceuyri$true.local.ancestry)),"\n")

## Three ancestries
## Not run: 
data(ceuchdyri)
res.eila <- eila(admixed  = ceuchdyri$admixed,
                 position = ceuchdyri$position,
                 anc1     = ceuchdyri$anc1,
                 anc2     = ceuchdyri$anc2,
                 anc3     = ceuchdyri$anc3)
cat("Overall accuracy:", mean(res.eila$local.ancestry ==
                                ceuchdyri$true.local.ancestry),"\n")

## End(Not run)