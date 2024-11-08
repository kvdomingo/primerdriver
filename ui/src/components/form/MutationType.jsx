import PropTypes from "prop-types";

function MutationType(props) {
  return (
    <div className="form-group">
      <label htmlFor="mutation_type">Mutation type</label>
      <select
        className="browser-default custom-select"
        id="mutation_type"
        name="mutation_type"
        onChange={e => props.handleChange(e.target.value)}
        required
        value={props.mutation_type}
      >
        <option value="">Select mutation type</option>
        <option value="sub">Substitution</option>
        <option value="ins">Insertion</option>
        <option value="del">Deletion</option>
      </select>
    </div>
  );
}

MutationType.propTypes = {
  handleChange: PropTypes.func,
  mutation_type: PropTypes.string,
};

export { MutationType };
