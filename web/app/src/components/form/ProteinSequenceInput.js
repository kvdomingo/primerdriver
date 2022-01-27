import PropTypes from "prop-types";
import { MDBIcon as Icon, MDBTooltip as Tooltip, MDBRow as Row, MDBCol as Col } from "mdbreact";

function ProteinSequenceInput(props) {
  return (
    <div className="form-group">
      <label htmlFor="sequence">Enter sequence</label>
      <Tooltip domElement tag="span" placement="right">
        <span>
          <Icon fas icon="question-circle" className="ml-2" />
        </span>
        <span>
          Simply type in or paste your primer amino acid sequence. Characters will automatically be filtered to show
          only standard one-letter amino acid IUPAC codes.
        </span>
      </Tooltip>
      <textarea
        className="form-control md-textarea"
        id="sequence"
        name="sequence"
        value={props.sequence}
        onChange={props.handleChangeSequence}
        onKeyUp={props.handleChangeSequence}
        onMouseUp={props.handleChangeSequence}
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

ProteinSequenceInput.propTypes = {
  cursorPosition: PropTypes.number.isRequired,
  handleChangeSequence: PropTypes.func.isRequired,
  sequence: PropTypes.string.isRequired,
  sequenceLength: PropTypes.number.isRequired,
};

export { ProteinSequenceInput };
