/*
 * Copyright (C) 2006 
 * Thomas van Dijk
 * Jan-Pieter van den Heuvel
 * Wouter Slob
 * 
 * Changed to TreeDecomposer by Matthieu Barba 2015
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
 
package nl.uu.cs.treewidth;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

import nl.uu.cs.treewidth.algorithm.Constructive;
import nl.uu.cs.treewidth.algorithm.GreedyDegree;
import nl.uu.cs.treewidth.algorithm.GreedyFillIn;
import nl.uu.cs.treewidth.algorithm.LexBFS;
import nl.uu.cs.treewidth.algorithm.MaximumCardinalitySearch;
import nl.uu.cs.treewidth.algorithm.MaximumCardinalitySearchMinimal;
import nl.uu.cs.treewidth.algorithm.PermutationGuesser;
import nl.uu.cs.treewidth.algorithm.PermutationToTreeDecomposition;
import nl.uu.cs.treewidth.algorithm.UpperBound;
import nl.uu.cs.treewidth.graph.Graph;
import nl.uu.cs.treewidth.graph.TDBag;
import nl.uu.cs.treewidth.input.DgfReader;
import nl.uu.cs.treewidth.input.GraphInput;
import nl.uu.cs.treewidth.input.InputException;
import nl.uu.cs.treewidth.input.GraphInput.InputData;
import nl.uu.cs.treewidth.ngraph.NGraph;
import nl.uu.cs.treewidth.ngraph.NTDBag;
import nl.uu.cs.treewidth.output.DotWriter;
import nl.uu.cs.treewidth.output.NeatoViewer;

/**
 * Program to get a tree decomposition of a graph
 */
public class TreeDecomposer {

	/**
	 * @param args 
	 */
	public static void main(String[] args) {
		// Get algos
		Constructive<InputData> theAlgorithm = null;
		Constructive[] algos = {
			new PermutationToTreeDecomposition<InputData>( new GreedyDegree<InputData>() ),
			new PermutationToTreeDecomposition<InputData>( new GreedyFillIn<InputData>() ),
			new PermutationToTreeDecomposition<InputData>( new LexBFS<InputData>() ),
			new PermutationToTreeDecomposition<InputData>( new MaximumCardinalitySearch<InputData>() ),
			//new PermutationToTreeDecomposition<InputData>( new MaximumCardinalitySearchMinimal<InputData>() ),
			new PermutationGuesser()
		};

		
		// Get exactly 3 arguments
		if( args.length < 3 ) {
			System.out.println("Need 3 args: method number, input_file.dgf, output_file.dot");
			System.out.println("Available methods:");
			// Print algos!
			int i = 1;
			for( Constructive a : algos ) {
				System.out.println( i++ + ") " + a.getName() );
			}
			return;
		}
		String method = args[0];
		String graphFile = args[1];
		String treeFile = args[2];
		String graphDotFile = args.length > 3 ? args[3] : "";
		
		// First control all three arguments
		// TODO
		
		// Then import the graph
		GraphInput in = new DgfReader( graphFile );
		NGraph<InputData> g;
		try {
			g = in.get();
		} catch (InputException e) {
			System.out.println( "Error opening file" );
			e.printStackTrace();
			return;
		}
		
		//Output comments
		String comments = g.getComments();
		if(comments.length()>0) {
			System.out.println( "The graph comments:");
			System.out.println( g.getComments() );
			System.out.println( "\nEnd of graph comments.");
		}
		
		int selected = 0;
		try {
			selected = Integer.parseInt(method)-1;
		} catch (Exception e){
			//Check input
			String error = "* ERROR: \"" + method + "\" is not a valid input *";
			String stars = "";
			for(int er=0;er<error.length();er++)
				stars = stars + "*";
			
			System.err.println("\r\n"+stars);
			System.err.println(error);
			System.err.println(stars+"\r\n");
			try {
				//Some delay to prevent the next loop from starting too soon.
				Thread.sleep(100);
			} catch (InterruptedException e1) {
				e1.printStackTrace();
			}
			main(args);
			return;
		}
		
		if( selected>=algos.length ) {
			System.out.println("Not in the list, fool." );
			System.exit(1);
		}
		
		// Selected algorithm
		theAlgorithm = algos[selected];

		// Load graph and run
		theAlgorithm.setInput(g);
		theAlgorithm.run();
		// Retrieve the Tree Decomposition
		NGraph<NTDBag<InputData>> td = theAlgorithm.getDecomposition();
		
		// Convert the Tree Decomposition to dot format
		String dotG = DotWriter.format( g );
		String dotTD = DotWriter.formatTD( td );
		// Write tree to file
		try {
			PrintWriter out = new PrintWriter( treeFile );
			out.println( dotTD );
			out.close();
		} catch(IOException e) {
			e.printStackTrace();
			return;	
		}
		
		// Write graph to file
		if (graphDotFile != "") {
			try {
				PrintWriter out = new PrintWriter( graphDotFile );
				out.println( dotG );
				out.close();
			} catch(IOException e) {
				e.printStackTrace();
				return;	
			}
		}
		/*
		NeatoViewer.present( dotG, graphFile+" : Input", 0, 0, false, true );
		
		String dotTD = DotWriter.formatTD( td );
		String title = graphFile+" : Treedecomposition by " + theAlgorithm.getName();
		if( theAlgorithm instanceof UpperBound ) {
			title = title + ", width " + ((UpperBound)theAlgorithm).getUpperBound();
		}
		NeatoViewer.present( dotTD, title , 500,0, false, true );
		*/
		
	}

}
