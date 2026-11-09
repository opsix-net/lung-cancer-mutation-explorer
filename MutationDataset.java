import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.TreeSet;

/** Load, validate, search, and compare the frozen lung cancer dataset. */
public class MutationDataset {
    public static final List<String> GENES = Collections.unmodifiableList(
        Arrays.asList("BRAF", "EGFR", "KEAP1", "KRAS", "STK11", "TP53"));
    private final Set<String> samples = new TreeSet<>();
    private final List<Mutation> mutations = new ArrayList<>();
    private final Set<String> types = new TreeSet<>();

    /** Stop on bad input data so incomplete files cannot silently bias counts. */
    public MutationDataset(Path directory) throws IOException {
        List<String> lines = Files.readAllLines(
            directory.resolve("samples.tsv"), StandardCharsets.UTF_8);
        if (lines.isEmpty() || !lines.get(0).equals("sample_id")) {
            throw new IOException("samples.tsv must begin with sample_id.");
        }
        for (int i = 1; i < lines.size(); i++) {
            String id = lines.get(i);
            if (!id.matches("[A-Za-z0-9-]+") || !samples.add(id)) {
                throw new IOException("Invalid or duplicate sample at line " + (i + 1));
            }
        }
        if (samples.isEmpty()) {
            throw new IOException("The sample list is empty.");
        }

        lines = Files.readAllLines(directory.resolve("mutations.tsv"), StandardCharsets.UTF_8);
        String header = "gene\tsample_id\tmutation_type\tprotein_change\tchromosome"
            + "\tstart\tend\treference\talternate";
        if (lines.isEmpty() || !lines.get(0).equals(header)) {
            throw new IOException("Unexpected mutations.tsv header.");
        }
        Set<String> identities = new HashSet<>();
        for (int i = 1; i < lines.size(); i++) {
            String[] fields = lines.get(i).split("\t", -1);
            try {
                if (fields.length != 9) {
                    throw new IllegalArgumentException("Expected nine columns.");
                }
                for (String field : fields) {
                    if (field.isBlank() || !field.equals(field.trim())) {
                        throw new IllegalArgumentException("Empty or padded field.");
                    }
                }
                if (!GENES.contains(fields[0]) || !samples.contains(fields[1])) {
                    throw new IllegalArgumentException("Unknown gene or sample.");
                }
                Mutation mutation = new Mutation(fields);
                if (mutation.start < 1 || mutation.end < mutation.start) {
                    throw new IllegalArgumentException("Invalid genomic coordinates.");
                }
                if (!identities.add(mutation.identity())) {
                    throw new IllegalArgumentException("Duplicate variant in a sample.");
                }
                mutations.add(mutation);
                types.add(mutation.type);
            } catch (IllegalArgumentException error) {
                throw new IOException("Invalid mutation at line " + (i + 1)
                    + ": " + error.getMessage(), error);
            }
        }
        // Sorting once gives stable results across runs and input row order.
        mutations.sort(Comparator.comparing((Mutation m) -> m.gene)
            .thenComparing(m -> m.sampleId).thenComparing(m -> m.chromosome)
            .thenComparingInt(m -> m.start).thenComparingInt(m -> m.end)
            .thenComparing(m -> m.reference).thenComparing(m -> m.alternate));
    }

    /** Accept gene symbols without requiring a particular capitalization. */
    public static String normalizeGene(String input) {
        String gene = input.trim().toUpperCase(Locale.ROOT);
        if (!GENES.contains(gene)) {
            throw new IllegalArgumentException("Choose a bundled gene: " + String.join(", ", GENES));
        }
        return gene;
    }

    /** Match a known type, or accept ALL to include every recorded type. */
    public String normalizeType(String input) {
        String value = input.trim();
        if (value.equalsIgnoreCase("ALL")) {
            return "ALL";
        }
        for (String type : types) {
            if (type.equalsIgnoreCase(value)) {
                return type;
            }
        }
        throw new IllegalArgumentException("Choose a listed mutation type or ALL.");
    }

    /** Return all matching records; the display limit is applied afterwards. */
    public List<Mutation> search(String gene, String type) {
        gene = normalizeGene(gene);
        type = normalizeType(type);
        List<Mutation> matches = new ArrayList<>();
        for (Mutation mutation : mutations) {
            if (mutation.gene.equals(gene)
                    && (type.equals("ALL") || mutation.type.equals(type))) {
                matches.add(mutation);
            }
        }
        return matches;
    }

    /** A sample can contain multiple mutations and must count only once. */
    public static Set<String> sampleIds(List<Mutation> records) {
        Set<String> ids = new HashSet<>();
        for (Mutation mutation : records) {
            ids.add(mutation.sampleId);
        }
        return ids;
    }

    /** Return first-only, second-only, both, and neither-recorded sample counts. */
    public int[] compare(String first, String second) {
        first = normalizeGene(first);
        second = normalizeGene(second);
        if (first.equals(second)) {
            throw new IllegalArgumentException("Choose two different genes.");
        }
        Set<String> a = sampleIds(search(first, "ALL"));
        Set<String> b = sampleIds(search(second, "ALL"));
        int[] counts = new int[4];
        for (String sample : samples) {
            if (a.contains(sample) && b.contains(sample)) {
                counts[2]++;
            } else if (a.contains(sample)) {
                counts[0]++;
            } else if (b.contains(sample)) {
                counts[1]++;
            } else {
                counts[3]++;
            }
        }
        return counts;
    }

    public int sampleCount() { return samples.size(); }
    public int mutationCount() { return mutations.size(); }
    public Set<String> mutationTypes() { return Collections.unmodifiableSet(types); }
}
