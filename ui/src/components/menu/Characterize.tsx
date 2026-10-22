import { useNavigate } from "@tanstack/react-router";
import { MDBCol as Col, MDBRow as Row } from "mdbreact";
import { type FormEvent, type SyntheticEvent, useEffect, useState } from "react";
import api from "../../api";
import { type ResultData, usePrimerDriverStore } from "../../store";
import { DnaSequenceInput, Form, MutationType, NumberMismatch } from "../form";
import LoadingScreen from "../shared/LoadingScreen";

function Characterize() {
  const [formData, setFormData] = useState({});
  const [sequence, setSequence] = useState("");
  const [mutationType, setMutationType] = useState("");
  const [mismatchedBases, setMismatchedBases] = useState(0);
  const [cursorPosition, setCursorPosition] = useState(1);
  const [sequenceLength, setSequenceLength] = useState(0);
  const [isValid, setIsValid] = useState(false);
  const [loading, setLoading] = useState(false);
  const navigate = useNavigate();
  const { setMode, setLoadedResultsFromModule, setResults } = usePrimerDriverStore();
  const mode = "CHAR" as const;

  useEffect(() => {
    setMode(mode);
  }, [setMode, mode]);

  function handleChangeSequence(e: SyntheticEvent<HTMLTextAreaElement>) {
    const sequence = e.currentTarget.value.toUpperCase().split("");
    const filteredSequence = [];

    for (const char of sequence) {
      if (["A", "G", "T", "C"].includes(char)) filteredSequence.push(char);
    }

    setSequence(filteredSequence.join(""));
    setCursorPosition(e.currentTarget.selectionStart + 1);
    setSequenceLength(e.currentTarget.value.length);
  }

  function handleSubmit(e: FormEvent) {
    e.preventDefault();
    setLoading(true);
    let data: ResultData;
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
        setLoadedResultsFromModule(true);
        setResults({ data, loaded: true });
        navigate({ to: "/results" });
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

export default Characterize;
