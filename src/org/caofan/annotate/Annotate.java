package org.caofan.annotate;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class Annotate {

	static final int splicing_region = 2;
	static final int nearthres = 1000;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		
		CommandLine cmd = parseArgs(args);
		
		long sum = 0;
		
		for(int i = 0; i< 1; i++) {
			long time1 = System.currentTimeMillis();
			
			HashMap<String, Attr> attrs = LoadReference.loadAttrs(cmd.getOptionValue('a'));

			//printAttrs(attrs);
			HashMap<String, ArrayList<Record>> reference = LoadReference
					.loadReference(cmd.getOptionValue('r'));
			//printRef(reference);
			
			HashMap<String, ArrayList<Record>> reject = null;
			if ( cmd.hasOption('b') ) {
				reject = LoadReference.loadRejectRegions( cmd.getOptionValue('b') );
			}
			
			int rejectedCount = 0;

			//printRef(reference);
			System.out.println("Sorted");
			try {
				for (String f : cmd.getOptionValues('i')) {
					System.out.println("\n=====================================");
					System.out.println("Processing " + f);
					System.out.println("=====================================");
					QueryFileReader reader = new QueryFileReader(f, cmd.getOptionValue('f').toUpperCase());
					PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(f + ".annotation."+cmd.getOptionValue('f').toLowerCase())));
					out.println(reader.getHeader());
					System.out.println("create output");
					QueryRecord r = reader.next();
					int count = 0;
					while(r != null) {
						if ( !isReject (reject, r) ) {

							count += 1;
							annotate(reference, r, attrs);
							out.println(r);
						} else{
							rejectedCount += 1;
						}
						r = reader.next();
					}
					System.out.println("Total records processed: " + (count + rejectedCount));
					System.out.println("total records annotated: " + count );
					System.out.println("Total records rejected: " + rejectedCount);
					out.close();
					reader.close();
					
				}
				
			} catch (IOException e) {
				System.err.println(e.getMessage());
				System.exit(1);
			}
			long time2 = System.currentTimeMillis();
			System.out.println(time2 - time1);
			sum += time2 - time1;
		}
		System.out.println("Average: " + sum*1.0/1);
	}
	
	private static Options createOptions() {
		Options opts = new Options();
		opts.addOption( OptionBuilder.withLongOpt( "reference" )
									 .withDescription( "the Gencode reference that could be obtained from UCSC." )
									 .hasArg()
									 .withArgName( "reference" )
									 .isRequired()
									 .create('r') );
		
		opts.addOption( OptionBuilder.withLongOpt( "attrs" )
									 .withDescription( "the Gencode attrs file that corresponds to the reference. ")
									 .isRequired()
									 .hasArg()
									 .withArgName( "Attr" )
									 .create('a') );
		
		opts.addOption( OptionBuilder.withLongOpt( "inputs" )
									 .withDescription( "input files" )
									 .hasArgs()
									 .withArgName( "Inputs" )
									 .isRequired()
									 .create('i') );
		
		opts.addOption( OptionBuilder.withLongOpt( "filetype" )
									 .withDescription("the type of the input files. " +
									 		"It can be either BED, GFF, MACS, or MACSSUMMIT." +
									 		"MACSSUMMIT is to use the summit of the MACS output for annotation." +
									 		"MACS and MACSSUMMIT uses the xls file output by MACS.")
									 .isRequired()
									 .hasArg()
									 .withArgName( "type" )
									 .create( 'f' ) );
		
		opts.addOption( OptionBuilder.withLongOpt( "rejectRegions" )
									 .withDescription( "Peaks fall in these regions will be ignored. Should be in BED format." )
									 .hasArgs()
									 .withArgName( "rejected" )
									 .create('b'));
		return opts;
	}
	
	public static CommandLine parseArgs(String[] args) {
		Options opts = createOptions();
		CommandLineParser parser = new GnuParser();
		
		try {
			CommandLine line = parser.parse(opts, args);
			return line;
		}
		catch( ParseException exp) {
			System.err.println( "Unexpected exception:" + exp.getMessage() );
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("ant", opts);
			System.exit(1);
			return null;
		}
		
		
	}
	
	
	public static void printAttrs(HashMap<String, Attr> map) {
		Iterator it = map.entrySet().iterator();
		while (it.hasNext()) {
			Map.Entry pairs = (Map.Entry)it.next();
			//System.out.println(pairs.getKey());
			Attr a = (Attr) pairs.getValue();
			System.out.println(a);
		}
	}
	

	
	
	public static void printRef(HashMap<String, ArrayList<Reference>> ref) {
		Iterator it = ref.entrySet().iterator();
		while (it.hasNext()) {
			Map.Entry pairs = (Map.Entry)it.next();
			//System.out.println(pairs.getKey());
			ArrayList<Reference> a = (ArrayList<Reference>) pairs.getValue();
			for(int i = 0; i < a.size(); i++) {
				System.out.println(a.get(i));
			}
		}
	}
	
	public static boolean isReject( HashMap< String, ArrayList< Record> > rejectRegion, QueryRecord q) {
		if (rejectRegion == null ) {
			return false;
		}
		ArrayList< Record > cRej = rejectRegion.get( q.getChrom() );
		if ( cRej == null ) {
			return false;
		} else {
			for ( Record r : cRej ) {
				if( overlapS( r.getStart(), r.getEnd(), q.getStart(), q.getEnd()) ) {
					return true;
				}
			}
		}
		
		return false;
	}
	
	public static void annotate(HashMap<String, ArrayList<Record>> ref, QueryRecord q, HashMap<String, Attr> attrs) {
		ArrayList<Record> cRef = ref.get(q.getChrom());
		if (cRef == null) {
			q.setAnn("None");
			q.setGenes("None");
		} else {
			HashSet<GeneAndType> exonic = new HashSet<GeneAndType>();
			HashSet<GeneAndType> splicing = new HashSet<GeneAndType>();
			HashSet<GeneAndType> intronic = new HashSet<GeneAndType>();
			HashSet<GeneAndType> utr5 = new HashSet<GeneAndType>();
			HashSet<GeneAndType> utr3 = new HashSet<GeneAndType>();
			HashSet<GeneAndType> upstream = new HashSet<GeneAndType>();
			HashSet<GeneAndType> downstream = new HashSet<GeneAndType>();
			boolean coding[] = {false, false, false, false, false, false, false};
			//For intergenic annotations
			String genel="None", gener="None";
			int distl = Integer.MAX_VALUE, distr = Integer.MAX_VALUE; 
			boolean foundGenic = false;
			boolean dummy;
			for (Record t : cRef) {
				Reference tr = (Reference) t;
				GeneAndType gt = new GeneAndType(tr.getName(), tr.isCoding()); //Currently we are not getting the distances to tss/tes here.
				if(overlapS(t.getStart(), t.getEnd(), q.getStart(), q.getEnd())) {
					foundGenic = true;
					for (int i = 0; i < tr.getExons().size(); i++) {
						Exon c = tr.getExons().get(i);
						if (overlapS(c.getStart(), c.getEnd(), q.getStart(),
								q.getEnd())) {
							//Whether the feature is in UTR regions.
							if (q.getEnd() < tr.getCdsStart()) {
								dummy = t.getStrand() == '+'? utr5.add(gt) : utr3.add(gt);   //Store the transcript name to make sure uniqueness. Look up other information from the attr.
							} else if(q.getStart() >= tr.getCdsEnd()) {
								dummy = t.getStrand() == '+'? utr3.add(gt) : utr5.add(gt);
							} else{
								exonic.add(gt);
								break;
							}
						} else if((i == 0 && overlapS(c.getStart(), c.getEnd() + splicing_region, q.getStart(), q.getEnd())) //The splicing regions
								|| (i == tr.getExons().size() - 1 && overlapS(c.getStart() - splicing_region, c.getEnd(), q.getStart(), q.getEnd())) 
								|| (i > 0 && i < tr.getExons().size() -1 && overlapS(c.getStart() - splicing_region, c.getEnd() + splicing_region, q.getStart(), q.getEnd()))) {
							splicing.add(gt);
						} else if(i < tr.getExons().size() - 1 && q.getStart() > c.getEnd() && q.getEnd() < tr.getExons().get(i+1).getStart()) {
							intronic.add(gt);
						}
					}
				} else if(!foundGenic && overlapS(t.getStart() - nearthres, t.getEnd(), q.getStart(), q.getEnd())) {
					dummy = t.getStrand() == '+' ? upstream.add(gt) : downstream.add(gt);
				} else if(!foundGenic && overlapS(t.getStart(), t.getEnd() + nearthres, q.getStart(), q.getEnd())) {
					dummy = t.getStrand() == '+' ? downstream.add(gt) : upstream.add(gt);
				} else if(!foundGenic) {
					if (t.getStart() >= q.getEnd() && t.getStart() - q.getEnd() + 1 < distr) {
						distr = t.getStart() - q.getEnd() + 1;
						gener = tr.getName();
					} else if (t.getEnd() <= q.getStart() && q.getStart() - t.getEnd() + 1 < distl) {
						distl = q.getStart() - t.getEnd() + 1;
						genel = tr.getName();
					}
				}
			}
			if(exonic.size() > 0) {
				q.genes = createOutput(exonic, attrs);
				q.ann = "Exonic";
			} else if(splicing.size() > 0) {
				q.genes = createOutput(splicing, attrs);
				q.ann = "Splicing";
			} else if(utr5.size() > 0 || utr3.size() > 0) {
				createOutputTwo(utr5, utr3, "UTR5", "UTR3",  attrs, q);
			} else if (intronic.size() > 0 ) {
				q.genes = createOutput(intronic, attrs);
				q.ann = "Intronic";
			} else if (upstream.size() > 0 || downstream.size() > 0) {
				createOutputTwo(upstream, downstream, "Upstream", "Downstream", attrs, q);
			} else {
				q.ann = "Intergenic";
				StringBuilder sb = new StringBuilder();
				if(!genel.equalsIgnoreCase("None")) {
					Attr assAttr = attrs.get(genel);
					sb.append(createIgStr(genel, distl, assAttr));
				} else {
					sb.append(genel);
				}
				sb.append(":");
				if(!gener.equalsIgnoreCase("None")) {
					Attr assAttr = attrs.get(genel);
					sb.append(createIgStr(gener, distr, assAttr));
				} else {
					sb.append(gener);
				}
				q.genes = sb.toString();
				
			}
		}
		
	}
	
	public static String createIgStr(String gene, int dist, Attr attr) {
		StringBuilder sb = new StringBuilder();
		if(attr != null) {
			sb.append(attr.geneName);
			sb.append(" (");
			sb.append(attr.geneType);
			sb.append(",");
		} else{
			sb.append(gene);
			sb.append(" (");
		}
		sb.append("dist=");
		sb.append(dist);
		sb.append(")");
		return sb.toString();
	}
		
	public static String createOutput(HashSet<GeneAndType> ann, HashMap<String, Attr> attrs) {
		boolean coding = false;
		for(GeneAndType a: ann) {
			if (a.isCoding) {
				coding = true;
				break;
			}
		}
		StringBuilder sb = new StringBuilder();
		HashSet<String> outStr = new HashSet<String>();
		for(GeneAndType a: ann) {
			if ((coding && a.isCoding) || !coding){
				Attr at = attrs.get(a.gene);
				sb = new StringBuilder();
				if (at != null) {
					sb.append(at.geneName);
					sb.append(" (");
					sb.append(at.geneType);
					sb.append(");");
				} else {
					sb.append(a.gene);
				}
				outStr.add(sb.toString());
			}
		}
		sb = new StringBuilder();
		for (String o: outStr) {
			sb.append(o);
		}
		return sb.toString();
	}
	
	public static void createOutputTwo(HashSet<GeneAndType> ann1, HashSet<GeneAndType> ann2, String type1, String type2, HashMap<String, Attr> attrs, QueryRecord q) {
		StringBuilder sb = new StringBuilder();
		if (ann1.size() > 0 && ann2.size() == 0) {
			q.ann =type1;
			q.genes = createOutput(ann1, attrs);
		} else if(ann1.size() == 0 && ann2.size() > 0) {
			q.ann = type2;
			q.genes = createOutput(ann2, attrs);
		} else if(ann1.size() > 0 && ann2.size() > 0) {
			q.ann = type1 + ":" + type2;
			sb.append(createOutput(ann1, attrs));
			sb.append(": ");
			sb.append(createOutput(ann2, attrs));
			q.genes = sb.toString();
		}
	}
	
	public static OverlapRegion overlap(int x1, int y1, int x2, int y2) {
		int minY = Math.min(y1, y2);
		int maxX = Math.max(x1, x2);
		return new OverlapRegion(maxX, minY);
	}
	
	public static boolean overlapS(int x1, int y1, int x2, int y2) {
		return Math.min(y1, y2) > Math.max(x1, x2);
	}
	
	static class GeneAndType {
		final String gene; //This should be a unique identifier of transcripts.
		final boolean isCoding; 
		int distTSS; //Distance to transcription starting site
		int distTES; //Distance to transcription ending site
		
		public GeneAndType(String g, boolean ic) {
			this(g, ic, Integer.MAX_VALUE, Integer.MAX_VALUE);
		}
		
		public GeneAndType(String g, boolean ic, int distTSS, int distTES) {
			this.gene = g;
			this.isCoding = ic;
			this.distTSS = distTSS;
			this.distTES = distTES;
		}
		
		public void setdTSS(int dist) {
			this.distTSS = dist;
		}
		
		public void setdTES(int dist) {
			this.distTES = dist;
		}
		
		@Override
		public boolean equals(Object obj) {
			if (obj == null) { return false; }
			if (obj == this) { return true; }
			if (obj.getClass() != getClass()) {
				return false;
			}
			
			GeneAndType g = (GeneAndType) obj;
			if (g.gene.equalsIgnoreCase(g.gene)) return true;
			return false;
		}
	}
	
	static class OverlapRegion {
		final int x;
		final int y;
		public OverlapRegion(int x, int y) {
			this.x = x;
			this.y = y;
		}
	}

}


