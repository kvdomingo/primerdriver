import { MDBCol as Col, MDBRow as Row } from "mdbreact";
import PropTypes from "prop-types";
import { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import api from "../../api/index.js";
import { usePrimerDriverContext } from "../../contexts/PrimerDriverContext.jsx";
import { DnaSequenceInput, Form, MutationType, NumberMismatch } from "../form";
import LoadingScreen from "../shared/LoadingScreen.jsx";

function Characterize(props) {
  const [formData, setFormData] = useState({});
  const [sequence, setSequence] = useState("");
  const [mutationType, setMutationType] = useState("");
  const [mismatchedBases, setMismatchedBases] = useState(0);
  const [cursorPosition, setCursorPosition] = useState(1);
  const [sequenceLength, setSequenceLength] = useState(0);
  const [isValid, setIsValid] = useState(false);
  const [loading, setLoading] = useState(false);
  const navigate = useNavigate();
  const { PDDispatch } = usePrimerDriverContext();
  const mode = "CHAR";

  useEffect(() => {
    PDDispatch({
      type: "updateMode",
      payload: mode,
    });
  }, [PDDispatch]);

  function handleChangeSequence(e) {
    const sequence = e.target.value.toUpperCase().split("");
    const filteredSequence = [];

    for (const char of sequence) {
      if (["A", "G", "T", "C"].includes(char)) filteredSequence.push(char);
    }

    setSequence(filteredSequence.join(""));
    setCursorPosition(e.target.selectionStart + 1);
    setSequenceLength(e.target.value.length);
  }

  function handleSubmit(e) {
    e.preventDefault();
    setLoading(true);
    let data;
    api.data
      .primerDriver(formData)
      .then(res => {
        data = res.data;
      })
      .catch(err => {
        console.log(err.message);
        data = "Request failed. Please try again later.";
      })
      .finally(() => {
        PDDispatch({
          type: "updateLoadedResults",
          payload: true,
        });
        PDDispatch({
          type: "updateResults",
          payload: {
            data,
            loaded: true,
          },
        });
        navigate("/results");
      });
  }

  function validateForm() {
    const validSequence = sequence.length > 0;
    const validMismatch = mismatchedBases > 0;
    const validMutation = mutationType !== "";
    setIsValid(validSequence && validMismatch && validMutation);
    setFormData({
      mode,
      sequence,
      mismatched_bases: mismatchedBases,
      mutation_type: mutationType,
    });
  }

  if (loading) return <LoadingScreen />;
  return (
    <Form
      title="Primer Characterization"
      handleValidate={validateForm}
      handleSubmit={handleSubmit}
      handleReset={props.handleReset}
      isValid={isValid}
    >
      <DnaSequenceInput
        cursorPosition={cursorPosition}
        handleChange={handleChangeSequence}
        sequence={sequence}
        sequenceLength={sequenceLength}
      />
      <Row className="row-cols-1 row-cols-md-2">
        <Col>
          <NumberMismatch
            handleChangeInt={setMismatchedBases}
            mismatched_bases={mismatchedBases}
          />
        </Col>
        <Col>
          <MutationType handleChange={setMutationType} mutation_type={mutationType} />
        </Col>
      </Row>
    </Form>
  );
}

Characterize.propTypes = {
  handleReset: PropTypes.func,
};

export default Characterize;
