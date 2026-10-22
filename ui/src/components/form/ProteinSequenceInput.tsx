import {
  MDBCol as Col,
  MDBIcon as Icon,
  MDBRow as Row,
  MDBTooltip as Tooltip,
} from "mdbreact";
import { type SyntheticEvent, useState } from "react";

interface ProteinSequenceInputProps {
  sequence: string;
  onChange: (value: string) => void;
}

function ProteinSequenceInput(props: ProteinSequenceInputProps) {
  const [cursorPosition, setCursorPosition] = useState(1);

  function handleChange(e: SyntheticEvent<HTMLTextAreaElement>) {
    props.onChange(e.currentTarget.value);
    setCursorPosition(e.currentTarget.selectionStart + 1);
  }

  return (
    <div className="form-group">
      <div className="alert alert-warning text-left" role="alert">
        <strong>Codon-optimized genes only.</strong> Protein-based design rebuilds a
        template by translating your sequence into the selected expression system&apos;s
        preferred codons, then designs against that reconstruction. The primers will
        anneal to your actual template only if the gene was synthesized with exactly
        those codons. For a cloned or cDNA-derived gene, use DNA-based design with the
        real coding sequence.
      </div>
      <label htmlFor="sequence">Enter sequence</label>
      <Tooltip domElement tag="span" placement="right">
        <span>
          <Icon fas icon="question-circle" className="ml-2" />
        </span>
        <span>
          Simply type in or paste your primer amino acid sequence. Characters will
          automatically be filtered to show only standard one-letter amino acid IUPAC
          codes.
        </span>
      </Tooltip>
      <textarea
        className="form-control md-textarea"
        id="sequence"
        name="sequence"
        value={props.sequence}
        onChange={handleChange}
        onKeyUp={handleChange}
        onMouseUp={handleChange}
        autoFocus
        required
        style={{ height: 128 }}
      />
      <Row className="row-cols-1 row-cols-md-2">
        <Col className="text-left">
          <small className="blue-grey-text">Cursor position: {cursorPosition}</small>
        </Col>
        <Col className="text-left text-md-right">
          <small className="blue-grey-text">
            Sequence length: {props.sequence.length}
          </small>
        </Col>
      </Row>
    </div>
  );
}

export { ProteinSequenceInput };
