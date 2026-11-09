/*
 * Lung Cancer Mutation Explorer. Java 17 or newer; no extra libraries.
 * In Binder: File > New > Terminal, then:
 *   mkdir -p /tmp/lung-explorer-classes
 *   javac -d /tmp/lung-explorer-classes *.java
 *   java -cp /tmp/lung-explorer-classes Main
 * Run from the project folder so the bundled data/ directory can be found.
 * Data: TCGA PanCancer Atlas via cBioPortal; see data/NOTICE.txt.
 */
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.TreeMap;

/** Console interface for exploring recorded mutations in a research cohort. */
public class Main {
    /** Start the application and report data failures without a stack trace. */
    public static void main(String[] args) {
        if (args.length > 1) {
            System.err.println("Usage: Main [data-directory]");
            System.exit(1);
        }
        try {
            MutationDataset data = new MutationDataset(
                Path.of(args.length == 1 ? args[0] : "data"));
            run(data, new Scanner(System.in));
        } catch (IOException error) {
            System.err.println("Cannot load dataset: " + error.getMessage());
            System.exit(1);
        }
    }

    /** Keep the menu available until the reader chooses Exit or input ends. */
    private static void run(MutationDataset data, Scanner input) {
        System.out.println("LUNG CANCER MUTATION EXPLORER");
        System.out.println("TCGA Lung Adenocarcinoma, PanCancer Atlas (2018)");
        System.out.println("Frozen six-gene snapshot. Research exploration only.");
        System.out.println("Samples: " + data.sampleCount()
            + "; mutation records: " + data.mutationCount());
        System.out.println("Counts describe samples, not patients or population risk.");
        System.out.println("No mutation recorded does not prove absence.");
        try {
            while (true) {
                System.out.println("\n1. Gene overview\n2. Search a gene"
                    + "\n3. Compare two genes\n4. Exit");
                int choice = readNumber(input, "Choice: ", 1, 4);
                if (choice == 4) {
                    break;
                } else if (choice == 1) {
                    overview(data);
                } else if (choice == 2) {
                    String gene = readGene(input, "Gene: ");
                    System.out.println("Types: ALL, " + String.join(", ", data.mutationTypes()));
                    String type;
                    while (true) {
                        System.out.print("Mutation type: ");
                        try {
                            type = data.normalizeType(input.nextLine());
                            break;
                        } catch (IllegalArgumentException error) {
                            System.out.println(error.getMessage());
                        }
                    }
                    int limit = readNumber(input, "Rows to display (1-20): ", 1, 20);
                    showSearch(data, gene, type, limit);
                } else {
                    String first = readGene(input, "First gene: ");
                    String second;
                    while (true) {
                        second = readGene(input, "Second gene: ");
                        if (!first.equals(second)) {
                            break;
                        }
                        System.out.println("Choose two different genes.");
                    }
                    showComparison(data, first, second);
                }
            }
        } catch (NoSuchElementException error) {
            // Piped input or a closed terminal can end at any prompt.
            System.out.println("\nInput ended.");
        }
        System.out.println("Goodbye.");
    }

    /** Read whole lines to avoid leftover newline problems after numbers. */
    private static int readNumber(Scanner input, String prompt, int low, int high) {
        while (true) {
            System.out.print(prompt);
            try {
                int value = Integer.parseInt(input.nextLine().trim());
                if (value >= low && value <= high) {
                    return value;
                }
            } catch (NumberFormatException error) {
                // Invalid integers are retried using the same bounded prompt.
            }
            System.out.println("Enter a whole number from " + low + " to " + high + ".");
        }
    }

    /** Keep unsupported genes separate from a valid search with no results. */
    private static String readGene(Scanner input, String prompt) {
        System.out.println("Genes: " + String.join(", ", MutationDataset.GENES));
        while (true) {
            System.out.print(prompt);
            try {
                return MutationDataset.normalizeGene(input.nextLine());
            } catch (IllegalArgumentException error) {
                System.out.println(error.getMessage());
            }
        }
    }

    /** Summarize each bundled gene using separate record and sample counts. */
    private static void overview(MutationDataset data) {
        System.out.println("Gene | Mutation records | Distinct samples with records");
        for (String gene : MutationDataset.GENES) {
            List<Mutation> records = data.search(gene, "ALL");
            System.out.println(gene + " | " + records.size() + " | "
                + MutationDataset.sampleIds(records).size());
        }
    }

    /** Summarize the full result before showing a limited set of detail rows. */
    private static void showSearch(MutationDataset data, String gene, String type, int limit) {
        List<Mutation> matches = data.search(gene, type);
        System.out.println("\nSearch: " + gene + " / " + type);
        System.out.println("Matching mutation records: " + matches.size());
        System.out.println("Distinct samples with matching records: "
            + MutationDataset.sampleIds(matches).size());
        if (matches.isEmpty()) {
            System.out.println("No mutation records match these filters in this snapshot.");
            return;
        }
        Map<String, Integer> counts = new TreeMap<>();
        for (Mutation mutation : matches) {
            counts.put(mutation.type, counts.getOrDefault(mutation.type, 0) + 1);
        }
        System.out.println("Mutation type counts: " + counts);
        System.out.println("Showing " + Math.min(limit, matches.size()) + " of "
            + matches.size() + " records. Coordinates use the source genome build GRCh37.");
        System.out.println("Sample | Type | Protein change | Genomic variant");
        for (int i = 0; i < Math.min(limit, matches.size()); i++) {
            System.out.println(matches.get(i));
        }
    }

    /** Partition the study sample list by recorded mutations in the two genes. */
    private static void showComparison(MutationDataset data, String first, String second) {
        int[] counts = data.compare(first, second);
        System.out.println("\nComparison: " + first + " and " + second + " (all recorded types)");
        System.out.println(first + " only: " + counts[0]);
        System.out.println(second + " only: " + counts[1]);
        System.out.println("Both genes: " + counts[2]);
        System.out.println("Neither gene has a mutation recorded: " + counts[3]);
        System.out.println("Total study samples: " + data.sampleCount());
        System.out.println("These are record-presence groups, not evidence of biological"
            + " absence or statistical association.");
    }
}
