import {
  MDBCol as Col,
  MDBIcon as Icon,
  MDBRow as Row,
  MDBTooltip as Tooltip,
} from "mdbreact";
import { type SyntheticEvent, useState } from "react";

interface DnaSequenceInputProps {
  sequence: string;
  onChange: (value: string) => void;
}

function DnaSequenceInput(props: DnaSequenceInputProps) {
  const [cursorPosition, setCursorPosition] = useState(1);

  function handleChange(e: SyntheticEvent<HTMLTextAreaElement>) {
    props.onChange(e.currentTarget.value);
    setCursorPosition(e.currentTarget.selectionStart + 1);
  }

  return (
    <div className="form-group">
      <label htmlFor="sequence">Enter sequence</label>
      <Tooltip domElement tag="span" placement="right">
        <span>
          <Icon fas icon="question-circle" className="ml-2" />
        </span>
        <span>
          Simply type in or paste your primer DNA sequence. Characters will
          automatically be filtered to show only A, T, C, G bases, and capitalized, as
          per IUPAC standards.
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

export { DnaSequenceInput };
