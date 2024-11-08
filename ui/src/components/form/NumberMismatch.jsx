import PropTypes from "prop-types";

function NumberMismatch(props) {
  return (
    <div className="form-group">
      <label htmlFor="mismatched_bases">Number of mismatched bases</label>
      <input
        className="form-control"
        min={0}
        type="number"
        id="mismatched_bases"
        name="mismatched_bases"
        value={props.mismatched_bases}
        onChange={e => props.handleChangeInt(Number.parseInt(e.target.value))}
        required
      />
    </div>
  );
}

NumberMismatch.propTypes = {
  handleChangeInt: PropTypes.func,
  mismatched_bases: PropTypes.number,
};

export { NumberMismatch };
