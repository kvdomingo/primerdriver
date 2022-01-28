import { useEffect, useState } from "react";
import { Image } from "cloudinary-react";
import { MDBContainer as Container } from "mdbreact";
import api from "../utils/Endpoints";

function Landing() {
  const [programVersion, setProgramVersion] = useState("");
  const [webVersion, setWebVersion] = useState("");
  const copyYear = new Date().getFullYear();

  useEffect(() => {
    api.data
      .version()
      .then(res => {
        setProgramVersion(res.data["program_version"]);
        setWebVersion(res.data["web_version"]);
      })
      .catch(err => console.log(err.message));
  }, []);

  return (
    <Container className="py-4 mx-0 px-4">
      <Image
        className="img-fluid text-center"
        cloudName="kdphotography-assets"
        publicId="primerdriver/PrimerDriver_logo"
        secure
        responsive
        responsiveUseBreakpoints
        dpr="auto"
        width="auto"
        crop="scale"
      />
      <div className="pt-5 text-center">
        <p>
          <b>PrimerDriver</b> ties together key primer design tools and protocols to automate the design of mutagenic
          PCR primers. This allows input of sequences from different sources. The tool can accommodate both DNA &
          protein sequences to incorporate base pair insertions, deletions, and substitutions as specified by the user.
        </p>
        <p>
          <b>PrimerDriver</b> can design primer pairs and compute for all oligonucleotide sequences that incorporate the
          desired mutations. This can cater to an array of primer designs from random mutations, site-directed single
          mutagenesis, batch design site-directed mutagenesis, to multiple-site mutagenesis.
        </p>
        <div className="text-center py-3">
          <a
            href="https://github.com/kvdomingo/primerdriver/releases"
            target="_blank"
            rel="noopener noreferrer"
            className="px-2"
          >
            Download CLI
          </a>
          {" | "}
          <a
            href="https://kvdomingo.github.io/primerdriver/"
            target="_blank"
            rel="noopener noreferrer"
            className="px-2"
          >
            Documentation
          </a>
        </div>
        <div className="text-center py-3 px-3">
          PrimerDriver v{programVersion} ({webVersion})
          <br />
          &copy; {copyYear} <a href="mailto:hello@kvdomingo.xyz">Kenneth Domingo</a> &amp;{" "}
          <a href="mailto:ngutierrez@evc.pshs.edu.ph">Nomer Gutierrez</a>
        </div>
      </div>
    </Container>
  );
}

export default Landing;
