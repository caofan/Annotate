package org.caofan.annotate;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;


public class QueryFileReader {
	private BufferedReader reader;
	private String filetype;
	private boolean reachEnd;
	private String header;
	public QueryFileReader(String filename, String filetype) {
		try {
			File infile = new File(filename);
			reader = new BufferedReader(new FileReader(infile.getAbsolutePath()));
			this.filetype = filetype;
			reachEnd = false;
			createHeader();
		} catch(IOException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}
	}
	
	private void createHeader() {
		StringBuilder sb = new StringBuilder();
		
		if(filetype.equalsIgnoreCase("macs") || filetype.equalsIgnoreCase("macsSummit")) {
			sb.append("chrom\t")
			  .append("start\t")
			  .append("end\t")
			  .append("length\t")
			  .append("strand\t")
			  .append("abs_summit\t")
			  .append("pileup\t")
			  .append("-log10(pvalue)\t")
			  .append("fold_enrichment\t")
			  .append("-log10(qvalue)\t")
			  .append("name\t")
			  .append("tagCount\t")
			  .append("RPK");
		} 
		else if(filetype.equalsIgnoreCase("gff")) {
			sb.append("seqname\t")
			  .append("source\t")
			  .append("feature\t")
			  .append("start\t")
			  .append("end\t")
			  .append("-log10(qvalue)\t")
			  .append("strand\t")
			  .append("frame\t")
			  .append("attribute");
		} 
		else if(filetype.equalsIgnoreCase("bed")) {
			sb.append("chrom\t")
			  .append("start\t")
			  .append("end\t")
			  .append("name\t")
			  .append("-log10(qvalue)\t")
			  .append("strand");
		}
		sb.append( "\ttype" );
		sb.append("\tgenes");
		header = sb.toString();
	}
	
	public String getHeader() {
		return header;
	}
	
	public void close() {
		try {
			reader.close();
		} catch(IOException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}
	}
	
	public QueryRecord next() {
		if (!reachEnd) {
			try {
				String line = reader.readLine();
				if(line == null) {
					reachEnd = true;
					return null;
				}
				else if(line.trim().equals("")) {
					return next();
				}
				
				//When passing the line, guarantee that it is valid.
				if(filetype.equalsIgnoreCase("gff")) {
					return fromGFF(line.trim());
				} 
				else if(filetype.equalsIgnoreCase("macs")) {
					return fromMACS(line.trim());
				} 
				else if(filetype.equalsIgnoreCase("bed")) {
					return fromBED(line.trim());
				}
				else if(filetype.equalsIgnoreCase("macsSummit")) {
					return fromMACSUseSummit(line.trim());
				}
				else {
					return fromOther(line.trim());
				}
			} catch(IOException e) {
				System.err.println(e.getMessage());
				return null;
			} 
		} else {
			return null;
		}
	}
	
	public QueryRecord fromGFF(String line) {
		StringTokenizer st = new StringTokenizer(line, "\t");
		String chrom = st.nextToken();
		for(int i = 0; i < 2; i++) {
			st.nextToken();
		}
		int start = Integer.parseInt(st.nextToken()) - 1;  //Convert to 0-based to be compatible with the reference.
		int end = Integer.parseInt(st.nextToken());
		double score = Double.parseDouble(st.nextToken());
		char strand = st.nextToken().charAt(0);
		return new QueryRecord(chrom, start, end, strand, score, line);
	}
	
	public QueryRecord fromMACS(String line) {
		//This method handles the xls files output by MACS
		if (line.startsWith("#") || line.equals("")) {
			return next();
		}
		StringTokenizer st = new StringTokenizer(line, "\t");
		String chrom = st.nextToken();
		String temp = st.nextToken();
		
		if ( ! isNumeric(temp)) {
			return next(); //Probably a header line
		}
		int start = Integer.parseInt(temp) - 1; //Convert to 0-based.
		int end = Integer.parseInt(st.nextToken());
		st.nextToken(); //The length
		char strand = st.nextToken().charAt(0);
		st.nextToken(); //The summit position
		double score = Double.parseDouble(st.nextToken());
		return new QueryRecord(chrom, start, end, strand, score, line);
	}
	
	public QueryRecord fromMACSUseSummit(String line) {
		// This method handles the xls files output by MACS
		// Returns the summits instead of the peaks.
		if (line.startsWith("#") || line.equals("")) {
			return next();
		}
		StringTokenizer st = new StringTokenizer(line, "\t");
		String chrom = st.nextToken();
		String temp = st.nextToken();

		if (!isNumeric(temp)) {
			return next(); // Probably a header line
		}
		st.nextToken();  //Skip the end of peak.
		st.nextToken(); // The length
		char strand = st.nextToken().charAt(0);
		temp = st.nextToken();
		String[] poses = temp.split(",");
		int start = Integer.parseInt(poses[0]) - 1; // The summit position, convert to 0-based
		int end = Integer.parseInt( poses[ poses.length - 1 ]) + 1;
		double score = Double.parseDouble(st.nextToken());
		return new QueryRecord(chrom, start, end, strand, score, line);
	}
	
	public QueryRecord fromBED(String line) {
		//Requires at least 3 fields for the BED format.
		if (line.startsWith("track") || line.startsWith("#")) {
			return next();
		}
		StringTokenizer st = new StringTokenizer(line, "\t");
		String chrom = st.nextToken();
		int start = Integer.parseInt(st.nextToken());
		int end = Integer.parseInt(st.nextToken());
		String name = null;
		double score = 0;
		char strand = '.';
		if (st.hasMoreTokens()) 
			name = st.nextToken();
		if (st.hasMoreTokens())
			score = Double.parseDouble(st.nextToken());
		if (st.hasMoreTokens())
			strand = st.nextToken().charAt(0);
		return new QueryRecord(chrom, start, end, strand, score, line);
	}
	
	public QueryRecord fromOther(String line) {
		//Process other file formats.
		//The coordinates are assumed to be 1-based inclusive.
		//Only require the first three columns as 
		//chrom\tstart\tend
		//header line should start with '#'
		if (line.startsWith("#")) {
			return next();
		}
		StringTokenizer st = new StringTokenizer(line,"\t");
		String  chrom = st.nextToken();
		int start = Integer.parseInt( st.nextToken() ) - 1;
		int end = Integer.parseInt( st.nextToken() );
		double score = 0;
		char strand = '.';
		return new QueryRecord(chrom, start, end, strand, score, line);
	}
	
	public static boolean isNumeric(String str) {
		try {
			int a = Integer.parseInt(str);
			return true;
		} catch( NumberFormatException nfe) {
			return false;
		}
	}
	
	
}
