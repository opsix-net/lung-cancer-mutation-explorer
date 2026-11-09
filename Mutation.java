/** One recorded genomic variant in one study sample. */
public class Mutation {
    public final String gene;
    public final String sampleId;
    public final String type;
    public final String protein;
    public final String chromosome;
    public final int start;
    public final int end;
    public final String reference;
    public final String alternate;

    /** Store an already validated row from the bundled tab-separated table. */
    public Mutation(String[] fields) {
        gene = fields[0];
        sampleId = fields[1];
        type = fields[2];
        protein = fields[3];
        chromosome = fields[4];
        start = Integer.parseInt(fields[5]);
        end = Integer.parseInt(fields[6]);
        reference = fields[7];
        alternate = fields[8];
    }

    /** Identify a variant within a sample without counting annotations twice. */
    public String identity() {
        return gene + "|" + sampleId + "|" + chromosome + "|" + start
            + "|" + end + "|" + reference + "|" + alternate;
    }

    /** Keep output compact while retaining the genomic coordinates and alleles. */
    public String toString() {
        return sampleId + " | " + type + " | " + protein + " | "
            + chromosome + ":" + start + "-" + end + " "
            + reference + ">" + alternate;
    }
}
