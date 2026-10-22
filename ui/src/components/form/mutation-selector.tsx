import { MDBCol as Col, MDBRow as Row } from "mdbreact";

interface MutationSelectorProps {
  handleChangeTarget: (value: string) => void;
  handleChangeReplacement: (value: string) => void;
  handleChangePosition: (value: number) => void;
  mutation_type: string;
  position: number;
  replacement: string;
  target: string;
}

function MutationSelector(props: MutationSelectorProps) {
  const columns = props.mutation_type === "sub" ? 3 : 2;

  if (props.mutation_type === "") return null;

  return (
    <Row className={`row-cols-1 row-cols-md-${columns}`}>
      {props.mutation_type !== "ins" ? (
        <Col className="form-group">
          <label htmlFor="target">Target</label>
          <input
            type="text"
            className="form-control"
            name="target"
            id="target"
            value={props.target}
            onChange={e => props.handleChangeTarget(e.target.value)}
            required
          />
        </Col>
      ) : null}
      <Col className="form-group">
        <label htmlFor="target">Position</label>
        <input
          type="number"
          min="0"
          className="form-control"
          name="position"
          id="position"
          value={props.position}
          onChange={e => props.handleChangePosition(Number(e.target.value))}
          required
        />
      </Col>
      {props.mutation_type !== "del" ? (
        <Col className="form-group">
          <label htmlFor="target">Replacement</label>
          <input
            type="text"
            className="form-control"
            name="replacement"
            id="replacement"
            min="0"
            value={props.replacement}
            onChange={e => props.handleChangeReplacement(e.target.value)}
            required={true}
          />
        </Col>
      ) : null}
    </Row>
  );
}

export { MutationSelector };
