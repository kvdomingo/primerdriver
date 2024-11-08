import {
  MDBCard as Card,
  MDBCardBody as CardBody,
  MDBCardHeader as CardHeader,
  MDBCol as Col,
  MDBCollapse as Collapse,
  MDBRow as Row,
} from "mdbreact";
import PropTypes from "prop-types";
import { useState } from "react";

function AdvancedSettings(props) {
  const [isOpen, setIsOpen] = useState(false);

  function toggleCollapse() {
    setIsOpen(!isOpen);
  }

  return (
    <div className="accordion my-3" id="accordionAdvanced">
      <Card className="bordered z-depth-0">
        <CardHeader id="headingAdvanced">
          <h5 className="mb-0 text-center">
            <button
              className="btn btn-link"
              color="none"
              id="show-advanced"
              type="button"
              data-toggle="collapse"
              data-target="#advanced"
              aria-expanded="true"
              aria-controls="advanced"
              onClick={toggleCollapse}
            >
              Advanced settings
            </button>
          </h5>
        </CardHeader>
        <Collapse
          id="advanced"
          aria-labelledby="headingAdvanced"
          data-parent="#accordionAdvanced"
          isOpen={isOpen}
        >
          <CardBody>
            <Row className="row-cols-3 text-center">
              <Col />
              <Col>Min</Col>
              <Col>Max</Col>
            </Row>
            <Row className="row-cols-3">
              <Col className="text-right">
                Melting point (<sup>o</sup>C)
              </Col>
              <Col>
                <div className="md-outline">
                  <input
                    type="number"
                    id="Tm_range_min"
                    name="Tm_range_min"
                    className="form-control"
                    value={props.tmRange[0]}
                    onChange={e => {
                      const { value } = e.target;

                      const newValue =
                        Number.parseFloat(value) < props.tmRange[1]
                          ? Number.parseFloat(value)
                          : props.tmRange[1] - 1;

                      props.setTmRange([newValue, props.tmRange[1]]);
                    }}
                  />
                </div>
              </Col>
              <Col>
                <div className="md-outline">
                  <input
                    type="number"
                    id="Tm_range_max"
                    name="Tm_range_max"
                    className="form-control"
                    value={props.tmRange[1]}
                    onChange={e => {
                      const { value } = e.target;

                      const newValue =
                        Number.parseFloat(value) > props.tmRange[0]
                          ? Number.parseFloat(value)
                          : props.tmRange[0] + 1;

                      props.setTmRange([props.tmRange[0], newValue]);
                    }}
                  />
                </div>
              </Col>
            </Row>
            <Row className="row-cols-3">
              <Col className="text-right">GC Content (%)</Col>
              <Col>
                <div className="md-outline">
                  <input
                    type="number"
                    id="gc_range_min"
                    name="gc_range_min"
                    className="form-control"
                    value={props.gcRange[0]}
                    onChange={e => {
                      const { value } = e.target;

                      const newValue =
                        Number.parseFloat(value) < props.gcRange[1]
                          ? Number.parseFloat(value)
                          : props.gcRange[1] - 1;

                      props.setGcRange([newValue, props.gcRange[1]]);
                    }}
                  />
                </div>
              </Col>
              <Col>
                <div className="md-outline">
                  <input
                    type="number"
                    id="gc_range_max"
                    name="gc_range_max"
                    className="form-control"
                    value={props.gcRange[1]}
                    onChange={e => {
                      const { value } = e.target;

                      const newValue =
                        Number.parseFloat(value) > props.gcRange[0]
                          ? Number.parseFloat(value)
                          : props.gcRange[0] + 1;

                      props.setGcRange([props.gcRange[0], newValue]);
                    }}
                  />
                </div>
              </Col>
            </Row>
            <Row className="row-cols-3">
              <Col className="text-right">Length (bp)</Col>
              <Col>
                <div className="md-outline">
                  <input
                    type="number"
                    id="length_min"
                    name="length_min"
                    className="form-control"
                    value={props.lengthRange[0]}
                    onChange={e => {
                      const { value } = e.target;

                      const newValue =
                        Number.parseFloat(value) < props.lengthRange[1]
                          ? Number.parseFloat(value)
                          : props.lengthRange[1] - 1;

                      props.setLengthRange([newValue, props.lengthRange[1]]);
                    }}
                  />
                </div>
              </Col>
              <Col>
                <div className="md-outline">
                  <input
                    type="number"
                    id="length_max"
                    name="length_max"
                    className="form-control"
                    value={props.lengthRange[1]}
                    onChange={e => {
                      const { value } = e.target;

                      const newValue =
                        Number.parseFloat(value) > props.lengthRange[0]
                          ? Number.parseFloat(value)
                          : props.lengthRange[0] + 1;

                      props.setLengthRange([props.lengthRange[0], newValue]);
                    }}
                  />
                </div>
              </Col>
            </Row>
            <Row className="row-cols-3">
              <Col className="text-right">5&apos; flanking region (bp)</Col>
              <Col>
                <div className="md-outline">
                  <input
                    type="number"
                    id="flank5_range_min"
                    name="flank5_range_min"
                    className="form-control"
                    value={props.flank5Range[0]}
                    onChange={e => {
                      const { value } = e.target;

                      const newValue =
                        Number.parseFloat(value) < props.flank5Range[1]
                          ? Number.parseFloat(value)
                          : props.flank5Range[1] - 1;

                      props.setFlank5Range([newValue, props.flank5Range[1]]);
                    }}
                  />
                </div>
              </Col>
              <Col>
                <div className="md-outline">
                  <input
                    type="number"
                    id="flank5_range_max"
                    name="flank5_range_max"
                    className="form-control"
                    value={props.flank5Range[1]}
                    onChange={e => {
                      const { value } = e.target;

                      const newValue =
                        Number.parseFloat(value) > props.flank5Range[0]
                          ? Number.parseFloat(value)
                          : props.flank5Range[0] + 1;

                      props.setFlank5Range([props.flank5Range[0], newValue]);
                    }}
                  />
                </div>
              </Col>
            </Row>
            <Row>
              <Col className="text-right">3&apos; flanking region (bp)</Col>
              <Col>
                <div className="md-outline">
                  <input
                    type="number"
                    id="flank3_range_min"
                    name="flank3_range_min"
                    className="form-control"
                    value={props.flank3Range[0]}
                    onChange={e => {
                      const { value } = e.target;

                      const newValue =
                        Number.parseFloat(value) < props.flank3Range[1]
                          ? Number.parseFloat(value)
                          : props.flank3Range[1] - 1;

                      props.setFlank3Range([newValue, props.flank3Range[1]]);
                    }}
                  />
                </div>
              </Col>
              <Col>
                <div className="md-outline">
                  <input
                    type="number"
                    id="flank3_range_max"
                    name="flank3_range_max"
                    className="form-control"
                    value={props.flank3Range[1]}
                    onChange={e => {
                      const { value } = e.target;

                      const newValue =
                        Number.parseFloat(value) > props.flank3Range[0]
                          ? Number.parseFloat(value)
                          : props.flank3Range[0] + 1;

                      props.setFlank3Range([props.flank3Range[0], newValue]);
                    }}
                  />
                </div>
              </Col>
            </Row>
            <hr />
            <Row className="row-cols-3">
              <Col className="text-right">5&apos; forward overlap (bp)</Col>
              <Col>
                <div className="md-outline">
                  <input
                    type="number"
                    id="forward_overlap5"
                    name="forward_overlap5"
                    className="form-control"
                    value={props.forwardOverlap5}
                    onChange={e =>
                      props.setForwardOverlap5(Number.parseInt(e.target.value))
                    }
                  />
                </div>
              </Col>
              <Col />
            </Row>
            <Row className="row-cols-3">
              <Col className="text-right">3&apos; forward overlap (bp)</Col>
              <Col>
                <div className="md-outline">
                  <input
                    type="number"
                    id="forward_overlap3"
                    name="forward_overlap3"
                    className="form-control"
                    value={props.forwardOverlap3}
                    onChange={e =>
                      props.setForwardOverlap3(Number.parseInt(e.target.value))
                    }
                  />
                </div>
              </Col>
              <Col />
            </Row>
            {(props.mode ?? "CHAR") === "PRO" ? (
              <Row className="row-cols-3">
                <Col className="text-right">Expression system</Col>
                <Col>
                  <select
                    className="browser-default custom-select"
                    id="expression_system"
                    name="expression_system"
                    value={props.expressionSystem}
                    onChange={e => props.setExpressionSystem(e.target.value)}
                  >
                    {(props.expressionList ?? []).map(exp => (
                      <option value={exp} key={exp}>
                        {exp}
                      </option>
                    ))}
                  </select>
                </Col>
              </Row>
            ) : null}
            <Row className="row-cols-3">
              <Col className="text-right">Primer type</Col>
              <Col>
                <select
                  className="browser-default custom-select"
                  id="primer_mode"
                  name="primer_mode"
                  value={props.primerMode}
                  onChange={e => props.setPrimerMode(e.target.value)}
                >
                  <option value="complementary">Complementary</option>
                  <option value="overlapping">Overlapping</option>
                </select>
              </Col>
              <Col />
            </Row>
            <Row className="row-cols-3">
              <Col />
              <Col>
                <div className="custom-control custom-checkbox mt-3 text-left">
                  <input
                    type="checkbox"
                    className="custom-control-input"
                    id="terminate_gc"
                    name="terminate_gc"
                    checked={props.terminateGc}
                    onChange={e => props.setTerminateGc(e.target.checked)}
                  />
                  <label className="custom-control-label" htmlFor="terminate_gc">
                    Terminates in G/C
                  </label>
                </div>
                <div className="custom-control custom-checkbox text-left">
                  <input
                    type="checkbox"
                    className="custom-control-input"
                    id="center_mutation"
                    name="center_mutation"
                    checked={props.centerMutation}
                    onChange={e => props.setCenterMutation(e.target.checked)}
                  />
                  <label className="custom-control-label" htmlFor="center_mutation">
                    Mutation at center of primer
                  </label>
                </div>
              </Col>
              <Col />
            </Row>
          </CardBody>
        </Collapse>
      </Card>
    </div>
  );
}

AdvancedSettings.propTypes = {
  tmRange: PropTypes.arrayOf(PropTypes.number),
  setTmRange: PropTypes.func,
  gcRange: PropTypes.arrayOf(PropTypes.number),
  setGcRange: PropTypes.func,
  lengthRange: PropTypes.arrayOf(PropTypes.number),
  setLengthRange: PropTypes.func,
  flank3Range: PropTypes.arrayOf(PropTypes.number),
  setFlank3Range: PropTypes.func,
  flank5Range: PropTypes.arrayOf(PropTypes.number),
  setFlank5Range: PropTypes.func,
  forwardOverlap3: PropTypes.number,
  setForwardOverlap3: PropTypes.func,
  forwardOverlap5: PropTypes.number,
  setForwardOverlap5: PropTypes.func,
  expressionSystem: PropTypes.string,
  setExpressionSystem: PropTypes.func,
  mode: PropTypes.string,
  primerMode: PropTypes.string,
  setPrimerMode: PropTypes.func,
  terminateGc: PropTypes.bool,
  setTerminateGc: PropTypes.func,
  centerMutation: PropTypes.bool,
  setCenterMutation: PropTypes.func,
  expressionList: PropTypes.array,
};

export { AdvancedSettings };
