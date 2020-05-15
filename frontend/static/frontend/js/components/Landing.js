import React, { Component } from 'react';
import { Image } from 'cloudinary-react';
import {
	MDBContainer as Container,
} from 'mdbreact';


export default class Landing extends Component {
	render() {
		return (
			<Container>
				<Image
					className='img-fluid px-md-5 pl-md-5 pr-md-5'
					cloudName='kdphotography-assets'
					publicId='primerdriver/PrimerDriver_logo'
					secure
					responsive
					responsiveUseBreakpoints
					dpr='auto'
					width='auto'
					crop='scale'
					/>
			    <div className="my-4">
			        <p>
			            <b>PrimerDriver</b> ties together key primer design tools and protocols to automate the design of mutagenic PCR primers. This allows input of sequences from different origins (i.e., FASTA, FASTQ, and manual input). The tool can accommodate both DNA & protein sequences to incorporate base pair insertions, deletions, and substitutions as specified by the user. <b>PrimerDriver</b> can design primer pairs and compute for all oligonucleotide sequences that incorporate the desired mutations. This can cater to an array of primer designs from random mutations, site-directed single mutagenesis, batch design site-directed mutagenesis, to multiple-site mutagenesis. Lastly, you can choose to download the results in different file formats - HTML, FASTA, and PDF.
			        </p>
			    </div>
			</Container>
		);
	}
}
