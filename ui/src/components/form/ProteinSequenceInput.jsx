import {
  MDBCol as Col,
  MDBIcon as Icon,
  MDBRow as Row,
  MDBTooltip as Tooltip,
} from "mdbreact";
import PropTypes from "prop-types";

function ProteinSequenceInput(props) {
  return (
    <div className="form-group">
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
        onChange={props.handleChangeSequence}
        onKeyUp={props.handleChangeSequence}
        onMouseUp={props.handleChangeSequence}
        autoFocus
        required
        style={{ height: 128 }}
      />
      <Row className="row-cols-1 row-cols-md-2">
        <Col className="text-left">
          <small className="blue-grey-text">
            Cursor position: {props.cursorPosition}
          </small>
        </Col>
        <Col className="text-left text-md-right">
          <small className="blue-grey-text">
            Sequence length: {props.sequenceLength}
          </small>
        </Col>
      </Row>
    </div>
  );
}

ProteinSequenceInput.propTypes = {
  cursorPosition: PropTypes.number,
  handleChangeSequence: PropTypes.func,
  sequence: PropTypes.string,
  sequenceLength: PropTypes.number,
};

export { ProteinSequenceInput };
