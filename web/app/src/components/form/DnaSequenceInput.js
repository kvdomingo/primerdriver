import PropTypes from "prop-types";
import { MDBIcon as Icon, MDBTooltip as Tooltip, MDBRow as Row, MDBCol as Col } from "mdbreact";

function DnaSequenceInput(props) {
  return (
    <div className="form-group">
      <label htmlFor="sequence">Enter sequence</label>
      <Tooltip domElement tag="span" placement="right">
        <span>
          <Icon fas icon="question-circle" className="ml-2" />
        </span>
        <span>
          Simply type in or paste your primer DNA sequence. Characters will automatically be filtered to show only A, T,
          C, G bases, and capitalized, as per IUPAC standards.
        </span>
      </Tooltip>
      <textarea
        className="form-control md-textarea"
        id="sequence"
        name="sequence"
        value={props.sequence}
        onChange={props.handleChange}
        onKeyUp={props.handleChange}
        onMouseUp={props.handleChange}
        autoFocus
        required
        style={{ height: 128 }}
      />
      <Row className="row-cols-1 row-cols-md-2">
        <Col className="text-left">
          <small className="blue-grey-text">Cursor position: {props.cursorPosition}</small>
        </Col>
        <Col className="text-left text-md-right">
          <small className="blue-grey-text">Sequence length: {props.sequenceLength}</small>
        </Col>
      </Row>
    </div>
  );
}

DnaSequenceInput.propTypes = {
  cursorPosition: PropTypes.number.isRequired,
  handleChange: PropTypes.func.isRequired,
  sequence: PropTypes.string.isRequired,
  sequenceLength: PropTypes.number.isRequired,
};

export { DnaSequenceInput };
